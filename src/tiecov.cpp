#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <set>

#include "commons.h"
#include <gclib/GArgs.h>
#include <gclib/GStr.h>
#include <gclib/GVec.hh>
#include "GSam.h"
#include <libBigWig/bigWig.h>

#define VERSION "0.0.6"

const char* USAGE="TieCov v" VERSION "\n"
"==================\n"
"The TieCov utility can take the output file produced by TieBrush and generate the following auxiliary files:\n"
" 1. BedGraph file with the coverage data\n"
" 2. Junction BED file\n"
" 3. a heatmap BED that uses color intensity to represent the number of samples that contain each position\n"
"==================\n"
"\n"
" usage: tiecov [-s out.sample] [-c out.coverage] [-j out.junctions] [-W] input\n"
"\n"
" Input arguments (required): \n"
"  input\t\talignment file in SAM/BAM/CRAM format\n"
"       "
"\n"
" Optional arguments (at least one of -s/-c/-j must be specified):\n"
"  -h,--help\tShow this help message and exit\n"
"  --version\tShow program version and exit\n"
"  -s\t\tBedGraph file with an estimate of the number of samples\n"
"    \t\twhich contain alignments for each interval.\n"
"  -c\t\tBedGraph (or BedWig with '-W') file with coverage\n"
"    \t\tfor all mapped bases.\n"
"  -j\t\tBED file with coverage of all splice-junctions\n"
"    \t\tin the input file.\n"
"  -W\t\tsave coverage in BigWig format. Default output\n"
"    \t\tis in Bed format\n";

GStr covfname, jfname, infname, sfname;
FILE* coutf=NULL;
FILE* joutf=NULL;
FILE* soutf=NULL;

GStr covfname_bw, jfname_bw, sfname_bw;
bigWigFile_t *coutf_bw = NULL;
bigWigFile_t *joutf_bw = NULL;
bigWigFile_t *soutf_bw = NULL;

std::vector<std::string> sample_info; // holds data about samples from the header

bool verbose=false;
bool bigwig=false;
int juncCount=0;

struct CJunc {
	int start, end;
	char strand;
	uint64_t dupcount;
	CJunc(int vs=0, int ve=0, char vstrand='+', uint64_t dcount=1):
	  start(vs), end(ve), strand(vstrand), dupcount(dcount) { }

	bool operator==(const CJunc& a) {
		return (strand==a.strand && start==a.start && end==a.end);
	}
//	bool operator<(const CJunc& a) { // sort no strand
//		if (start==a.start) return (end<a.end);
//		else return (start<a.start);
//	}

//    bool operator<(const CJunc& a) { // sort by strand first
//        if (strand==a.strand){
//            if (start==a.start){
//                return (end<a.end);
//            }
//            else{
//                return (start<a.start);
//            }
//        }
//        else{
//            return strand<a.strand;
//        }
//    }

    bool operator<(const CJunc& a) { // sort by strand last
        if (start==a.start){
            if(end==a.end){
                return strand<a.strand;
            }
            else{
                return (end<a.end);
            }
        }
        else{
            return (start<a.start);
        }
    }

	void add(CJunc& j) {
       dupcount+=j.dupcount;
	}

	void write(FILE* f, const char* chr) {
		juncCount++;
		fprintf(f, "%s\t%d\t%d\tJUNC%08d\t%ld\t%c\n",
				chr, start-1, end, juncCount, (long)dupcount, strand);
	}
};

GArray<CJunc> junctions(64, true);

void addJunction(GSamRecord& r, int dupcount) {
	char strand = r.spliceStrand();
//	if (strand!='+' && strand!='-') return; // TODO: should we output .?
	for (int i=1;i<r.exons.Count();i++) {
		CJunc j(r.exons[i-1].end+1, r.exons[i].start-1, strand,
				dupcount);
		int ei;
		int r=junctions.AddIfNew(j, &ei);
		if (r==-1) {//existing junction, update
			junctions[ei].add(j);
		}
	}
}

void flushJuncs(FILE* f, const char* chr) {
    for (int i=0;i<junctions.Count();i++) {
    	junctions[i].write(f,chr);
    }
    junctions.Clear();
    junctions.setCapacity(128);
}

void processOptions(int argc, char* argv[]);

void addSamples(GSamRecord& r, std::vector<int>& cur_samples, std::vector<std::set<int>>& bvec, int b_start){ // same as addMean below - but here it computes the actual number of samples when used with an index
    bam1_t* in_rec=r.get_b();
    int pos=in_rec->core.pos; // 0-based
    b_start--; //to make it 0-based
    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int oplen=bam_cigar_oplen(cigar_full[c]);
        switch(opcode){
            case BAM_CINS: // no change in coverage and position
                break;
            case BAM_CDEL: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CREF_SKIP: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CSOFT_CLIP:
                break;
            case BAM_CMATCH: // base match - add coverage
                for(int i=0;i<oplen;i++) {
                    bvec[pos-b_start].insert(cur_samples.begin(),cur_samples.end()); // add sample ids to the running set
                    pos++;
                }
                break;
            default:
                GError("ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
        }
    }
}

void addMean(GSamRecord& r, int val, std::vector<std::pair<float,uint64_t>>& bvec, int b_start){ // for YX (number of samples) we are not interested in the sum but rather the average number of smaples that describe the position. Giving a heatmap
    bam1_t* in_rec=r.get_b();
    int pos=in_rec->core.pos; // 0-based
    b_start--; //to make it 0-based
    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int oplen=bam_cigar_oplen(cigar_full[c]);
        switch(opcode){
            case BAM_CINS: // no change in coverage and position
                break;
            case BAM_CDEL: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CREF_SKIP: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CSOFT_CLIP:
                break;
            case BAM_CMATCH: // base match - add coverage
                for(int i=0;i<oplen;i++) {
                    bvec[pos-b_start].first+=(val-bvec[pos-b_start].first)/bvec[pos-b_start].second; // dynamically compute average
                    bvec[pos-b_start].second++;
                    pos++;
                }
                break;
            default:
                GError("ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
        }
    }
}

void move2bsam(std::vector<std::set<int>>& bvec_idx,std::vector<std::pair<float,uint64_t>>& bvec){
    for(int i=0;i<bvec_idx.size();i++){
        bvec[i].second=bvec_idx[i].size();
    }
}

//b_start MUST be passed 1-based
void addCov(GSamRecord& r, int val, GVec<uint64_t>& bvec, int b_start) {
	bam1_t* in_rec=r.get_b();
    int pos=in_rec->core.pos; // 0-based
    b_start--; //to make it 0-based
    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int oplen=bam_cigar_oplen(cigar_full[c]);
        switch(opcode){
            case BAM_CINS: // no change in coverage and position
                break;
            case BAM_CDEL: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CREF_SKIP: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CSOFT_CLIP:
                break;
            case BAM_CMATCH: // base match - add coverage
                for(int i=0;i<oplen;i++) {
                    bvec[pos-b_start]+=val;
                    pos++;
                }
                break;
            default:
                GError("ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
        }
    }
}

//b_start MUST be passed 1-based
void flushCoverage(FILE* outf,sam_hdr_t* hdr, GVec<uint64_t>& bvec,  int tid, int b_start) {
  if (tid<0 || b_start<=0) return;
  int i=0;
  b_start--; //to make it 0-based;
  while (i<bvec.Count()) {
     uint64_t ival=bvec[i];
     int j=i+1;
     while (j<bvec.Count() && ival==bvec[j]) {
    	 j++;
     }
     if (ival!=0){
         fprintf(outf, "%s\t%d\t%d\t%ld\n", hdr->target_name[tid], b_start+i, b_start+j, (long)ival);
     }
     i=j;
  }
}

void flushCoverage(bigWigFile_t* outf,sam_hdr_t* hdr, GVec<uint64_t>& bvec,  int tid, int b_start) {
    if (tid<0 || b_start<=0) return;
    int i=0;
    b_start--; //to make it 0-based;
    bool first = true;
    while (i<bvec.Count()) {
        uint64_t ival=bvec[i];
        int j=i+1;
        while (j<bvec.Count() && ival==bvec[j]) {
            j++;
        }
        if (ival!=0){
            char *chromsUse[] = {hdr->target_name[tid]};
            uint32_t bw_start[] = {(uint32_t)b_start+i};
            uint32_t bw_end[] = {(uint32_t)b_start+j};
            float values[] = {(float)ival};
            if(first){
                if(bwAddIntervals(outf, chromsUse, bw_start, bw_end,values, 1)){
                    std::cerr<<"error bw"<<std::endl;
                    exit(-1);
                }
                first=false;
            }
            else{
                if(bwAppendIntervals(outf, bw_start, bw_end,values, 1)){
                    std::cerr<<"error bw"<<std::endl;
                    exit(-1);
                }
            }
        }
        i=j;
    }
}

void flushCoverage(FILE* outf,sam_hdr_t* hdr, std::vector<std::pair<float,uint64_t>>& bvec,  int tid, int b_start) {
    if (tid<0 || b_start<=0) return;
    int i=0;
    b_start--; //to make it 0-based;
    while (i<bvec.size()) {
        uint64_t ival=bvec[i].second;
        float hval = bvec[i].first;
        int j=i+1;
        while (j<bvec.size() && ival==bvec[j].second) {
            j++;
        }
        if (ival!=0)
            fprintf(outf, "%s\t%d\t%d\t%ld\t%f\n", hdr->target_name[tid], b_start+i, b_start+j, (long)ival,hval);
        i=j;
    }
}

void discretize(std::vector<std::pair<float,uint64_t>>& bvec1){
    for(auto& val : bvec1){
        //val.second = std::ceil(val.first);
        val.second = std::ceil(val.first);
        val.first = 0;
    }
}

void average_sample(std::vector<uint64_t>& bvec,float thresh){
    // iterate
    // find min and max of the range of values
    // find percentage by which values can be similar to group together
    // find difference between two points
    // if within difference - keep average and compare next observation with an average
    // if not in range - write all previous values with the average and start again
    uint64_t max_val = *std::max_element(std::begin(bvec), std::end(bvec));
    uint64_t min_val = *std::min_element(std::begin(bvec), std::end(bvec));

    for(auto& val : bvec){

    }
}

void normalize(std::vector<std::pair<float,uint64_t>>& bvec,float mint, float maxt, int num_samples){ // normalizes values to a specified range
    float denom = num_samples;
    float mult = (maxt-mint);

    for (auto& val : bvec){
        val.first = (val.second/denom)*mult+mint;
    }
}

void load_sample_list(std::vector<int>& lst,std::string& sl_fname,std::vector<std::string>& sample_info){
    std::ifstream sl_fp(sl_fname, std::ios_base::binary);
    std::string line;
    std::set<std::string> sample_lst_set;
    std::set<std::string>::iterator sl_it;
    while (sl_fp >> line){
        sample_lst_set.insert(line);
    }

    for(int i=0;i<sample_info.size();i++){ // find positions in the index to be extracted
        sl_it = sample_lst_set.find(sample_info[i]);
        if(sl_it!=sample_lst_set.end()){ // found
            lst.push_back(i);
        }
    }

    sl_fp.close();
}

// >------------------ main() start -----
int main(int argc, char *argv[])  {
    processOptions(argc, argv);
    //htsFile* hts_file=hts_open(infname.chars(), "r");
    //if (hts_file==NULL)
    //   GError("Error: could not open alignment file %s \n",infname.chars());
	GSamReader samreader(infname.chars(), SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);

    if (!covfname.is_empty()) {
       if (covfname=="-" || covfname=="stdout")
    	   coutf=stdout;
       else {
           if(!bigwig){
               if(std::strcmp(covfname.substr(covfname.length()-9,9).chars(),".bedgraph")!=0){ // if name does not end in .bedgraph
                   covfname.append(".bedgraph");
               }
               coutf=fopen(covfname.chars(), "w");
               if (coutf==NULL) GError("Error creating file %s\n",
                                       covfname.chars());
               fprintf(coutf, "track type=bedGraph\n");
           }
           else{ // initialize bigwig
               if(std::strcmp(covfname_bw.substr(covfname_bw.length()-7,7).chars(),".bigwig")!=0){ // if name does not end in .bigwig
                   covfname_bw.append(".bigwig");
               }
               if(bwInit(1<<17) != 0) {
                   fprintf(stderr, "Received an error in bwInit\n");
                   return 1;
               }
               coutf_bw = bwOpen((char*)covfname_bw.chars(), NULL, "w");
               if(!coutf_bw) {
                   fprintf(stderr, "An error occurred while opening example_output.bw for writing\n");
                   return 1;
               }
               //Allow up to 10 zoom levels, though fewer will be used in practice
               if(bwCreateHdr(coutf_bw, 10)){
                   std::cerr<<"error bw"<<std::endl;
                   exit(-1);
               }
               //Create the chromosome lists
               // need to get from htslib
               char *chroms[samreader.header()->n_targets];
               uint32_t lens[samreader.header()->n_targets];
               for(int bwi=0;bwi<samreader.header()->n_targets;bwi++){
                   chroms[bwi] = samreader.header()->target_name[bwi];
                   lens[bwi] = samreader.header()->target_len[bwi];
               }
               coutf_bw->cl = bwCreateChromList(chroms, lens, samreader.header()->n_targets);
               if(!coutf_bw->cl){
                   std::cerr<<"error bw"<<std::endl;
                   exit(-1);
               }

               //Write the header
               if(bwWriteHdr(coutf_bw)){
                   std::cerr<<"error bw"<<std::endl;
                   exit(-1);
               }
           }
       }
    }
    if (!jfname.is_empty()) {
        if(std::strcmp(jfname.substr(jfname.length()-4,4).chars(),".bed")!=0){ // if name does not end in .bed
            jfname.append(".bed");
        }
       joutf=fopen(jfname.chars(), "w");
       if (joutf==NULL) GError("Error creating file %s\n",
        		  jfname.chars());
       fprintf(joutf, "track name=junctions\n");
    }
    if (!sfname.is_empty()) {
        if(std::strcmp(sfname.substr(sfname.length()-9,9).chars(),".bedgraph")!=0){ // if name does not end in .bedgraph
            sfname.append(".bedgraph");
        }
        soutf=fopen(sfname.chars(), "w");
        if (soutf==NULL) GError("Error creating file %s\n",
                                sfname.chars());
        fprintf(soutf, "track type=bedGraph name=\"Sample Count Heatmap\" description=\"Sample Count Heatmap\" visibility=full graphType=\"heatmap\" color=200,100,0 altColor=0,100,200\n");
    }

    int prev_tid=-1;
    GVec<uint64_t> bcov(2048*1024);
    std::vector<std::pair<float,uint64_t>> bsam(2048*1024,{0,1}); // number of samples. 1st - current average; 2nd - total number of values
    std::vector<std::set<int>> bsam_idx(2048*1024,std::set<int>{}); // for indexed runs
    int b_end=0; //bundle start, end (1-based)
    int b_start=0; //1 based
    GSamRecord brec;
    while (samreader.next(brec)) {
        uint32_t dupcount=0;
        std::vector<int> cur_samples;
        int endpos=brec.end;
        if (brec.refId()!=prev_tid || (int)brec.start>b_end) {
            if (prev_tid>=0) {
              if (coutf)
                  flushCoverage(coutf,samreader.header(), bcov, prev_tid, b_start);
              if(coutf_bw)
                  flushCoverage(coutf_bw,samreader.header(), bcov, prev_tid, b_start);
              if (soutf) {
                  discretize(bsam);
                  normalize(bsam,0.1,1.5,sample_info.size());
                  flushCoverage(soutf,samreader.header(),bsam,prev_tid,b_start);
              }
              if (joutf)
                  flushJuncs(joutf, samreader.refName(prev_tid));
            }
            b_start=brec.start;
            b_end=endpos;
            if (coutf || coutf_bw) {
                bcov.setCount(0);
                bcov.setCount(b_end-b_start+1);
            }
            if (soutf) {
                bsam.clear();
                bsam.resize(b_end-b_start+1,{0,1});
                bsam_idx.clear();
                bsam_idx.resize(b_end-b_start+1,std::set<int>{});
            }
            prev_tid=brec.refId();
        } else { //extending current bundle
            if (b_end<endpos) {
                b_end=endpos;
                bcov.setCount(b_end-b_start+1, (int)0);
                if (soutf){
                    bsam.resize(b_end-b_start+1,{0,1});
                    bsam_idx.resize(b_end-b_start+1,std::set<int>{});
                }
            }
        }
        int accYC = 0;
        accYC = brec.tag_int("YC", 1);
        if(coutf || coutf_bw){
            addCov(brec, accYC, bcov, b_start);
        }
        if (joutf && brec.exons.Count()>1) {
            addJunction(brec, accYC);
        }

        if(soutf){
            addSamples(brec,cur_samples,bsam_idx,b_start);
            float accYX = 0;
            accYX = (float)brec.tag_int("YX", 1);
            addMean(brec, accYX, bsam, b_start);
        }
	} //while GSamRecord emitted
	if (coutf) {
       flushCoverage(coutf,samreader.header(), bcov, prev_tid, b_start);
       if (coutf!=stdout) fclose(coutf);
	}
	if (soutf) {
        discretize(bsam);
        normalize(bsam,0.1,1.5,sample_info.size());
        flushCoverage(soutf,samreader.header(),bsam,prev_tid,b_start);
        if (soutf!=stdout) fclose(soutf);
	}
	if (joutf) {
		flushJuncs(joutf, samreader.refName(prev_tid));
		fclose(joutf);
	}

	// same for BigWig
    if (coutf_bw) {
        flushCoverage(coutf_bw,samreader.header(), bcov, prev_tid, b_start);
        bwClose(coutf_bw);
        bwCleanup();
    }
    if (soutf_bw) {
        bwClose(soutf_bw);
        bwCleanup();
    }
    if (joutf_bw) {
        bwClose(joutf_bw);
        bwCleanup();
    }

}// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;verbose;version;DVWhc:s:j:");
    args.printError(USAGE, true);
    if (args.getOpt('h') || args.getOpt("help")) {
        GMessage(USAGE);
        exit(1);
    }

    if (args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }

    if ((args.getOpt('c') || args.getOpt('s') || args.getOpt('j'))==0){
        GMessage(USAGE);
        GMessage("\nError: at least one of -c/-j/-s arguments required!\n");
        exit(1);
    }

    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    bigwig=args.getOpt('W')!=NULL;

    if (verbose) {
        fprintf(stderr, "Running TieCov " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    covfname=args.getOpt('c');
    jfname=args.getOpt('j');
    sfname=args.getOpt('s');

    covfname_bw=args.getOpt('c');
    jfname_bw=args.getOpt('j');
    sfname_bw=args.getOpt('s');

    if (args.startNonOpt()==0) {
        GMessage(USAGE);
        GMessage("\nError: no input file provided!\n");
        exit(1);
    }
    
    infname=args.nextNonOpt();
}
