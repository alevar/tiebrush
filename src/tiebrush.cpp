#include <vector>
#include <cstring>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

#include "commons.h"
#include "GSam.h"
#include <gclib/GArgs.h>
#include "tmerge.h"
#include <gclib/GBitVec.h>

#define VERSION "0.0.7"

const char* USAGE = "TieBrush v" VERSION "\n"
                              "==================\n"
                              "Summarize and filter read alignments from multiple sequencing samples "
                              "(taken as sorted SAM/BAM/CRAM files). This utility aims to merge/collapse "
                              "\"duplicate\" read alignments across multiple sequencing samples (inputs), "
                              "adding custom SAM tags in order to keep track of the \"alignment multiplicity\" "
                              "count (how many times the same alignment is seen across all input data) and "
                              "\"sample count\" (how many samples show that same alignment).\n"
                              "==================\n"
                              "\n usage: tiebrush  [-h] -o OUTPUT [-L|-P|-E] [-S] [-M] [-N max_NH_value] "
                              "[-Q min_mapping_quality] [-F FLAGS] ...\n"
                              "\n"
                              " Input arguments:\n"
                              "  ...  \t\t\tinput alignment files can be provided as a space-delimited \n"
                              "       \t\t\tlist of filenames or as a text file containing a list of\n"
                              "       \t\t\tfilenames, one per line\n"
                              "\n"
                              " Required arguments:\n"
                              "  -o\t\t\tFile for BAM output\n"
                              "\n"
                              " Optional arguments:\n"
                              "  -h,--help\t\tShow this help message and exit\n"
                              "  --version\t\tShow the program version and exit\n"
                              "  -L,--full\t\tIf enabled, only reads with the same CIGAR\n"
                              "           \t\tand MD strings will be grouped and collapsed.\n"
                              "           \t\tBy default, TieBrush will consider the CIGAR\n"
                              "           \t\tstring only when grouping reads\n"
                              "           \t\tOnly one of -L, -P or -E options can be enabled\n"
                              "  -P,--clip\t\tIf enabled, reads will be grouped by clipped\n"
                              "           \t\tCIGAR string. In this mode 5S10M5S and 3S10M3S\n"
                              "           \t\tCIGAR strings will be grouped if the coordinates\n"
                              "           \t\tof the matching substring (10M) are the same\n"
                              "           \t\tbetween reads\n"
                              "  -E,--exon\t\tIf enabled, reads will be grouped if their exon\n"
                              "           \t\tboundaries are the same. This option discards\n"
                              "           \t\tany structural variants contained in mapped\n"
                              "           \t\tsubstrings of the read and only considers start\n"
                              "           \t\tand end coordinates of each non-splicing segment\n"
                              "           \t\tof the CIGAR string\n"
                              "  -S,--keep-supp\tIf enabled, supplementary alignments will be\n"
                              "                \tincluded in the collapsed groups of reads.\n"
                              "                \tBy default, TieBrush removes any mappings\n"
                              "                \tnot listed as primary (0x800).\n"
                              "  --keep-secondary\tIf enabled, secondary alignments will be\n"
                              "                   \tincluded in the collapsed groups of reads.\n"
                              "                   \tBy default, TieBrush removes any mappings\n"
                              "                   \tnot listed as primary (0x100). Note, that if enabled,\n"
                              "                   \teach secondary mapping will count as a separate read\n"
                              "  -M,--keep-unmap\tIf enabled, unmapped reads will be retained (uncollapsed)\n"
                              "                 \tin the output. By default, TieBrush removes any\n"
                              "                 \tunmapped reads\n"
                              "  -N\t\t\tMaximum NH score of the reads to retain\n"
                              "  -Q\t\t\tMinimum mapping quality of the reads to retain\n"
                              "  -F\t\t\tBits in SAM flag to use in read comparison. Only reads that\n"
							  "    \t\t\thave specified flags will be merged together (default: 0)\n"
							  "  -A,--collapse-same\tIf enabled, will collapse same read alignment\n"
							  "    \tduplicated just for pairing reasons (default: disabled)\n"
							  "  --store-frac\t\tIf enabled, will store fractional values for YC tag\n"
							  "    \tusing NH tag (1.0/NH) for non-TieBrush records (default: disabled)\n";

// 1. add mode to select representative alignment
// 2. add mode to select consensus sequence
// 3. add mode (current) - random read selection
// 4. replace options with argparse?
// 5. in the help indicate what options are default
// 6. fix PG/RG sample confusion

enum TMrgStrategy {
	tMrgStratCIGAR=0,  // same CIGAR (MD may differ)
	tMrgStratFull, // same CIGAR and MD
	tMrgStratClip,   // same CIGAR after clipping
	tMrgStratExon    // same exons
};

struct Options{
    int max_nh = MAX_INT;
    int min_qual = -1;
    bool keep_unmapped = true;
    bool keep_supplementary = false;
    bool keep_secondary = false;
    uint32_t flags = 0;
    bool collapse_same = false;
    bool store_frac = false;
} options;

TMrgStrategy mrgStrategy=tMrgStratCIGAR;
TInputFiles inRecords;

GStr outfname;
GSamWriter* outfile=NULL;

uint64_t inCounter=0;
uint64_t outCounter=0;

bool verbose=false;

struct GSegNode {
	GSeg seg;
	GSegNode* next;
	inline uint start() { return seg.start; }
	inline uint end() { return seg.end; }
	GSegNode(uint segstart=0, uint segend=0, GSegNode* nnext=NULL):seg(segstart, segend),
			next(nnext) {}
	GSegNode(GSeg& exon, GSegNode* nnext=NULL):seg(exon),
			next(nnext) {}
};

struct GSegList { //per sample per strand
  GSegNode* startNode;
  uint last_pos;
  int last_dist;
  GSegList():startNode(NULL),last_pos(0),last_dist(-1) { }

  ~GSegList() {
	  clear();
  }

  void reset() {
	  clear();
	  last_pos=0;
	  last_dist=-1;
	  //startNode=new GSegNode();
	  startNode=NULL;
  }

  void clear() { //delete all nodes
	  GSegNode* p=startNode;
	  GSegNode* next;
	  while (p) {
		  next=p->next;
		  delete p;
		  p=next;
	  }
	  startNode=NULL;
  }

  void clearTo(GSegNode* toNode) {
	  //clear every node up to and *including* toNode
	  GSegNode* p=startNode;
	  GSegNode* next;
	  while (p && p!=toNode) {
		  next=p->next;
		  delete p;
		  p=next;
	  }
	  if (p==NULL) GError("Error: clearTo did not find target node %d-%d!\n",
			  toNode->start(), toNode->end());
	  next=toNode->next;
	  delete toNode;
	  startNode=next;
  }

 void mergeRead(GSamRecord& r) {
	 if (startNode==NULL) {
		 startNode=new GSegNode(r.exons[0]);
		 GSegNode* cn=startNode;
		 for (int i=1;i<r.exons.Count();i++){
			 GSegNode* n=new GSegNode(r.exons[i]);
			 cn->next=n;
			 cn=n;
		 }
		 return;
	 }
	 GSegNode *n=startNode;
	 GSegNode *prev=NULL;
	 for (int i=0;i<r.exons.Count();i++) {
		 GSeg& e=r.exons[i];
		 while (n) {
            if (e.end < n->start()) {
              //exon should be inserted before n!
              GSegNode* nw=new GSegNode(e, n);
              if (n==startNode)
            	startNode=nw;
              else
                prev->next=nw;
              //n is unchanged
          	  prev=nw;
              break; //check next exon against n
            }
            // e.end >= node start
            if (e.start <= n->end()) { //overlap!
              //union should replace/change n
              if (e.start<n->start()) n->seg.start=e.start;
              if (e.end>n->end()) n->seg.end=e.end;
              //we have to check if next nodes are now overlapped as well
              GSegNode* next=n->next;
              while (next && next->start()<=n->end()) {
            	  //overlap, swallow this next node
            	  uint nend=next->end();
            	  n->next=next->next;
            	  delete next;
            	  next=n->next;
                  if (nend>n->end()) {
                	n->seg.end=nend;
                	break; //there cannot be overlap with next node
                  }
              } //while overlap with next nodes
              break; //check next exon against this extended n
            }
            // e.start > node end
            prev=n;
            n=n->next;
		 }
	 }
 }

 int processRead(GSamRecord& r) { //return d=current bundle extent upstream
	 //if the read starts after a gap, d=0
	 //this should only be called ONCE per collapsed read and sample
	 //should NOT be called on reads coming from already merged samples!
	 if (last_pos==r.start) { //already called on the same sample and start position
	     mergeRead(r);
		 return last_dist;
	 }
	 int d=0;
	 GSegNode* node=startNode;
	 GSegNode* prev=NULL;
	 while (node && node->start()< r.start) {
		 prev=node;
		 node=node->next;
	 }
	 //prev is the last segment starting before r
	 if (prev) {
		 if (prev->end()>=r.start)  // r overlaps prev segment
			d=r.start - prev->start();
		 if (d==0)
			clearTo(prev); //clear all nodes including prev
	 }

     if (last_pos!=r.start) {
    	 last_pos=r.start;
    	 last_dist=d;
     }
     mergeRead(r);
	 return d;
 }


};


struct RDistanceData {
  GVec<GSegList> fsegs; //forward strand segs for each sample
  GVec<GSegList> rsegs; //reverse strand segs for each sample
  void init(int num_samples) {
	  fsegs.Resize(num_samples);
	  rsegs.Resize(num_samples);
	  this->reset();
  }
  void reset() {
	  for (int i=0;i<fsegs.Count();i++) {
           fsegs[i].reset();
           rsegs[i].reset();
	  }
  }
};

RDistanceData rspacing;

// check the two reads for compatibility with user provided flags
int cmpFlags(GSamRecord& a, GSamRecord& b){
    if(options.flags == 0){
        return 0;
    }
    if((options.flags & a.get_b()->core.flag) == (options.flags & b.get_b()->core.flag)){ // make sure the user-requested flags are the same between reads
        return 0;
    }
    return 1;
}

int cmpFull(GSamRecord& a, GSamRecord& b) {
    //-- Flags
    bool cmp_flag = cmpFlags(a,b);
    if(cmp_flag!=0) return cmp_flag;
	//-- CIGAR && MD strings
	if (a.get_b()->core.n_cigar!=b.get_b()->core.n_cigar) return ((int)a.get_b()->core.n_cigar - (int)b.get_b()->core.n_cigar);
	int cigar_cmp=0;
	if (a.get_b()->core.n_cigar>0)
		cigar_cmp=memcmp(bam_get_cigar(a.get_b()) , bam_get_cigar(b.get_b()), a.get_b()->core.n_cigar*sizeof(uint32_t) );
	if (cigar_cmp!=0) return cigar_cmp;
	// compare MD tag
	char* aMD=a.tag_str("MD");
	char* bMD=b.tag_str("MD");
	if (aMD==NULL || bMD==NULL) {
		if (aMD==bMD) return 0;
		if (aMD!=NULL) return 1;
		return -1;
	}
    return strcmp(aMD, bMD);
}

int cmpCigar(GSamRecord& a, GSamRecord& b) {
    bool cmp_flag = cmpFlags(a,b);
    if(cmp_flag!=0) return cmp_flag;
	if (a.get_b()->core.n_cigar!=b.get_b()->core.n_cigar) return ((int)a.get_b()->core.n_cigar - (int)b.get_b()->core.n_cigar);
	if (a.get_b()->core.n_cigar==0) return 0;
	return memcmp(bam_get_cigar(a.get_b()) , bam_get_cigar(b.get_b()), a.get_b()->core.n_cigar*sizeof(uint32_t) );
}

int cmpCigarClip(GSamRecord& a, GSamRecord& b) {
    bool cmp_flag = cmpFlags(a,b);
    if(cmp_flag!=0) return cmp_flag;
	uint32_t a_clen=a.get_b()->core.n_cigar;
	uint32_t b_clen=b.get_b()->core.n_cigar;
	uint32_t* a_cstart=bam_get_cigar(a.get_b());
	uint32_t* b_cstart=bam_get_cigar(b.get_b());
	while (a_clen>0 &&
			((*a_cstart) & BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) { a_cstart++; a_clen--; }
	while (a_clen>0 &&
			(a_cstart[a_clen-1] & BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) a_clen--;
	while (b_clen>0 &&
			((*b_cstart) & BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) { b_cstart++; b_clen--; }
	while (b_clen>0 &&
			(b_cstart[b_clen-1] & BAM_CIGAR_MASK)==BAM_CSOFT_CLIP) b_clen--;
	if (a_clen!=b_clen) return (int)a_clen-(int)b_clen;
	if (a_clen==0) return 0;
	return memcmp(a_cstart, b_cstart, a_clen*sizeof(uint32_t));
}

int cmpExons(GSamRecord& a, GSamRecord& b) {
    bool cmp_flag = cmpFlags(a,b);
    if(cmp_flag!=0) return cmp_flag;
	if (a.exons.Count()!=b.exons.Count()) return (a.exons.Count()-b.exons.Count());
	for (int i=0;i<a.exons.Count();i++) {
		if (a.exons[i].start!=b.exons[i].start)
			return ((int)a.exons[i].start-(int)b.exons[i].start);
		if (a.exons[i].end!=b.exons[i].end)
			return ((int)a.exons[i].end-(int)b.exons[i].end);
	}
	return 0;
}


//keep track of all SAM alignments starting at the same coordinate
// that were merged into a single alignment
class SPData { // Same Point data
    bool settled; //real SP data, owns its r data and deallocates it on destroy
  public:
	double accYC; //if any merged records had YC tag, their values are accumulated here
	int64_t accYX; //if any merged records had YX tag, their values are accumulated here
	int64_t maxYD; //max bundle extent upstream in all samples (0 when no preceding overlapping reads)
	GBitVec* samples; //which samples were the collapsed ones coming from
	                 //number of bits set will be stored as YX:i:(samples.count()+accYX)
    std::vector<uint64_t> sample_dupcounts; // within each sample how many records were collapsed
                                            // number of bits set will be stored as YX:i:(samples.count()+accYX)
	int dupCount; //duplicity count - how many single-alignments were merged into r
	              // will be stored as tag YC:i:(dupCount+accYC)
	GSamRecord* r;
	char tstrand; //'-','+' or '.'
    SPData(GSamRecord* rec=NULL):settled(false), accYC(0.0), accYX(0), maxYD(0),samples(NULL),
    		dupCount(0), r(rec), tstrand('.') {
    	if (r!=NULL) tstrand=r->spliceStrand();
    }

    ~SPData() {
    	if (settled && r!=NULL) delete r;
    	if (samples!=NULL) delete samples;
        if (!sample_dupcounts.empty()) sample_dupcounts.clear();
    }

    void detach(bool dontFree=true) { settled=!dontFree; }
    //detach(true) must never be called before settle()

    void settle(TInputRecord& trec) { //becomes a standalone SPData record
    	// duplicates the current record
    	settled=true;
    	trec.disown();
    	if (samples==NULL) {
            samples = new GBitVec(inRecords.freaders.Count());
        }
        if (sample_dupcounts.empty()){
            sample_dupcounts.assign(inRecords.freaders.Count(),0);
        }

    	if (trec.tbMerged) {
    		// For TieBrush records, YC should already be a double value
    		accYC=r->tag_float("YC");
    		if (accYC == 0.0) accYC = 1.0; // fallback if tag not found
    		accYX=r->tag_int("YX", 1);
    		maxYD=r->tag_int("YD", 0);
    	} else {
    		if (options.store_frac) {
    			int nh = r->tag_int("NH", 1);
    			accYC = 1.0 / nh;
    		} else {
    			accYC = 1.0;
    		}
    		++dupCount;
    		samples->set(trec.fidx);
            sample_dupcounts[trec.fidx]++;
    	}
    }

    void dupAdd(TInputRecord& trec) { //merge an external SAM record into this one
    	if (!settled) GError("Error: cannot merge a duplicate into a non-settled SP record!\n");
    	GSamRecord& rec=*trec.brec;
    	//WARNING: rec MUST be a "duplicate" of current record r
    	if (trec.tbMerged) {
    		// For TieBrush records, YC should already be a double value
    		double yc_val = rec.tag_float("YC");
    		if (yc_val == 0.0) yc_val = 1.0; // fallback if tag not found
    		accYC += yc_val;
    		accYX+=rec.tag_int("YX", 1);
    		int64_t vYD=rec.tag_int("YD",0);
    		if (vYD>maxYD) maxYD=vYD; //keep only maximum YD value
    	} else {
    		//avoid collapsing same read alignment duplicated just for pairing reasons (optional)
    		if (!options.collapse_same || 
    		    !samples->test(trec.fidx) || rec.pairOrder()!=r->pairOrder() ||
    				strcmp(r->name(), rec.name())!=0) {
    		   if (options.store_frac) {
    			   int nh = rec.tag_int("NH", 1);
    			   accYC += 1.0 / nh;
    		   } else {
    			   accYC += 1.0;
    		   }
    		   dupCount++;
    		   samples->set(trec.fidx);
    		   sample_dupcounts[trec.fidx]++;
    		}
    	}
    }

    bool operator<(const SPData& b) {
    	if (r==NULL || b.r==NULL) GError("Error: cannot compare uninitialized SAM records\n");
    	if (r->refId()!=b.r->refId()) return (r->refId()<b.r->refId());
    	//NOTE: already assuming that start&end must match, no matter the merge strategy
    	if (r->start!=b.r->start) return (r->start<b.r->start);
    	if (tstrand!=b.tstrand) return (tstrand<b.tstrand);
    	if (r->end!=b.r->end) return (r->end<b.r->end);
    	if (mrgStrategy==tMrgStratFull) {
    		int ret=cmpFull(*r, *(b.r));
    		return (ret<0);
    	}
    	else
    		switch (mrgStrategy) {
    		  case tMrgStratCIGAR: return (cmpCigar(*r, *(b.r))<0); break;
    		  case tMrgStratClip: return (cmpCigarClip(*r, *(b.r))<0); break;
    		  case tMrgStratExon: return (cmpExons(*r, *(b.r))<0); break;
    		  default: GError("Error: unknown merge strategy!\n");
    		}
    	return false;
    }

    bool operator==(const SPData& b) {
    	if (r==NULL || b.r==NULL) GError("Error: cannot compare uninitialized SAM records\n");
    	if (r->refId()!=b.r->refId() || r->start!=b.r->start || tstrand!=b.tstrand ||
    			r->end!=b.r->end) return false;
    	if (mrgStrategy==tMrgStratFull) return (cmpFull(*r, *(b.r))==0);
    	else
    		switch (mrgStrategy) {
    		  case tMrgStratCIGAR: return (cmpCigar(*r, *(b.r))==0); break;
    		  case tMrgStratClip: return (cmpCigarClip(*r, *(b.r))==0); break;
    		  case tMrgStratExon: return (cmpExons(*r, *(b.r))==0); break;
    		  default: GError("Error: unknown merge strategy!\n");
    		}
    	return false;
    }
};

void processOptions(int argc, char* argv[]);

void addPData(TInputRecord& irec, GList<SPData>& spdlst) {
  //add and collapse if match found
	SPData* newspd=new SPData(irec.brec);
	if (spdlst.Count()>0) {
		//find if irec can merge into existing SPData
		SPData* spf=spdlst.AddIfNew(newspd, false);
		if (spf!=newspd) { //matches existing SP entry spf, merge
		    // TODO: consensus can be collected in dupAdd by remembering all MD data
		    //   initiate array of length (sequence) and store most common base
		    //   - use raw seq data in the record (already array)
		    //   - bit operations?
		    //     - write a separate function for this - that way can implement a simple approach first (plain iteration) and then improve it
			spf->dupAdd(irec); //update existing SP entry
			delete newspd;
			return;
		} // not a novel SP data
		//else newspd was added as a separate entry of spdlst
	}
	else { // empty list, just add this
		spdlst.Add(newspd);
	}
	newspd->settle(irec); //keep its own SAM record copy
}

void flushPData(GList<SPData>& spdlst){ //write spdata to outfile
  if (spdlst.Count()==0) return;
  // write SAM records in spdata to outfile
  for (int i=0;i<spdlst.Count();++i) {
	  SPData& spd=*(spdlst.Get(i));
	  int64_t accYX=spd.accYX;
	  int dSamples=spd.samples->count(); //this only has direct, non-TieBrush samples
	  accYX+=dSamples;
	  spd.r->add_double_tag("YC", spd.accYC);
	  spd.r->add_int_tag("YX", (accYX == 1) ? 1 : accYX);
	  int dmax=spd.maxYD;
	  for(int s=spd.samples->find_first();s>=0;s=spd.samples->find_next(s)) {
	    	if (spd.tstrand=='+' || spd.tstrand=='.') {
	    	   int r=rspacing.fsegs[s].processRead(*spd.r);
	    	   if (r>dmax) dmax=r;
	    	}
	    	if (spd.tstrand=='-' || spd.tstrand=='.') {
	    	   int r=rspacing.rsegs[s].processRead(*spd.r);
	    	   if (r>dmax) dmax=r;
	    	}
	  } //for each bit index/sample
	  spd.maxYD=dmax;
	  if (spd.maxYD>0) spd.r->add_int_tag("YD", spd.maxYD);
	  else spd.r->remove_tag("YD");
	  outfile->write(spd.r);

	  outCounter++;
  }
  spdlst.Clear();
}

bool passes_options(GSamRecord* brec){
    if(!options.keep_supplementary && brec->get_b()->core.flag & 0x100) return false;
    if(!options.keep_secondary && brec->get_b()->core.flag & 0x200) return false;
    if(!options.keep_unmapped && brec->isUnmapped()) return false;
    if(brec->mapq()<options.min_qual) return false;
    int nh = brec->tag_int("NH");
    if(nh>options.max_nh)  return false;

    return true;
}

// >------------------ main() start -----

// merging indices can be done as follows:
// 1. as we parse the alignments from each input - check duplicity for each sample
// or perhaps it can be a separate method
// 1. the next index can not have fewer entries - only more
// 2. so what we need to do is insert 0s where appropriate
// 3.    everything that is a 0 in a new file is a new 0
// 4.    everything that is >0 in a new file - is the next old value

// reserve space for a header - always constant, this way we can simply replace values there and not write anything else
// process indices one by one
//  1.

int main(int argc, char *argv[])  {
	inRecords.setup(VERSION, argc, argv);
	processOptions(argc, argv);
	int numSamples=inRecords.start();
	outfile=new GSamWriter(outfname, inRecords.header(), GSamFile_BAM);
	rspacing.init(numSamples);
	TInputRecord* irec=NULL;
	GSamRecord* brec=NULL;

	GList<SPData> spdata(true, true, true); //list of Same Position data, with all possibly merged records
	bool newChr=false;
	int prev_pos=-1;
	int prev_tid=-1;
	while ((irec=inRecords.next())!=NULL) {
		 brec=irec->brec;
         if(!passes_options(brec)) continue;
		 inCounter++;
		 int tid=brec->refId();
		 int pos=brec->start; //1-based

		 if (tid!=prev_tid) {
			 if (prev_tid!=-1) newChr=true;
			 prev_tid=tid;
			 prev_pos=-1;
		 }
		 if (pos!=prev_pos) { //new position
			 flushPData(spdata); //also adds read data to rspacing
			 prev_pos=pos;
		 }
		 if (newChr) {
			 rspacing.reset();
			 newChr=false;
		 }
		 addPData(*irec, spdata);
	}
    flushPData(spdata);
	inRecords.stop();

    delete outfile;

    //if (verbose) {
    double p=100.00 - (double)(outCounter*100.00)/(double)inCounter;
    GMessage("%ld input records written as %ld (%.2f%% reduction)\n", inCounter, outCounter, p);
    //}
}
// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;full;clip;exon;keep-supp;keep-secondary;keep-unmap;collapse-same;store-frac;SMLPEDVho:N:Q:F:A");
    args.printError(USAGE, true);

    if (args.getOpt('h') || args.getOpt("help")) {
        fprintf(stdout,"%s",USAGE);
        exit(0);
    }

    if (args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }

    if (args.startNonOpt()==0) {
        GMessage(USAGE);
        GMessage("\nError: no input provided!\n");
        exit(1);
    }
    outfname=args.getOpt('o');
    if (outfname.is_empty()) {
        GMessage(USAGE);
        GMessage("\nError: output filename must be provided (-o)!\n");
        exit(1);
    }

    GStr max_nh_str=args.getOpt('N');
    if (!max_nh_str.is_empty()) {
        options.max_nh=max_nh_str.asInt();
    }
    GStr min_qual_str=args.getOpt('Q');
    if (!min_qual_str.is_empty()) {
        options.min_qual=min_qual_str.asInt();
    }
    GStr flag_str=args.getOpt('F');
    if (!flag_str.is_empty()) {
        options.flags=flag_str.asInt();
    }
    options.keep_supplementary = (args.getOpt("keep-supp")!=NULL || args.getOpt("S")!=NULL);
    options.keep_secondary = (args.getOpt("keep-secondary")!=NULL);
    options.keep_unmapped = (args.getOpt("keep-unmap")!=NULL || args.getOpt("M")!=NULL);
    options.collapse_same = (args.getOpt("collapse-same")!=NULL || args.getOpt('A')!=NULL);
    options.store_frac = (args.getOpt("store-frac")!=NULL);

    bool stratF=(args.getOpt("full")!=NULL || args.getOpt('L')!=NULL);
    bool stratP=(args.getOpt("clip")!=NULL || args.getOpt('P')!=NULL);
    bool stratE=(args.getOpt("exon")!=NULL || args.getOpt('E')!=NULL);
    if (stratF | stratP | stratE) {
        if (!(stratF ^ stratP ^ stratE))
            GError("Error: only one merging strategy can be requested.\n");
        if (stratF) mrgStrategy=tMrgStratFull;
        else if (stratP) mrgStrategy=tMrgStratClip;
        else mrgStrategy=tMrgStratExon;
    }

    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    if (verbose) {
        fprintf(stderr, "Running TieBrush " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    const char* ifn=NULL;
    while ( (ifn=args.nextNonOpt())!=NULL) {
        //input alignment files
        std::string absolute_ifn = get_full_path(ifn);
        inRecords.addFile(absolute_ifn.c_str());
    }
}