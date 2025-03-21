#include <vector>
#include <limits>
#include <tuple>
#include <fstream>
#include <functional>
#include <numeric>
#include <cmath>
#include <map>
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


#define VERSION "0.0.3"

const char* USAGE="TieCov v" VERSION " usage:\n"
                  " tiecov [-b out.flt.bam] [-c out.bedgraph] [-j out.junctions.bed] in.bam\n"
                  " Other options: \n"
                  "  -N   : maximum NH score (if available) to include when reporting coverage\n"
                  "  -Q   : minimum mapping quality to include when reporting coverage\n"
                  "  -s   : compute coverage across splice junctions and write to file\n"
                  "  -l   : read length - if specified will be used in place of the max_start & max end as a hard cap\n";

struct Filters{
    int max_nh = MAX_INT;
    int min_qual = -1;
    int max_len = MAX_INT;
} filters;

struct OutputFMT{
    bool report_sjs = false;
    std::string sjs_fname  = "";
    std::string sjs_combined = "";
} outFmt;

GStr covfname, jfname, bfname, infname;
FILE* boutf=NULL;
FILE* coutf=NULL;
FILE* joutf=NULL;
FILE* sjs_strand_outf=NULL;
FILE* sjs_comb_outf=NULL;

bool debugMode=false;
bool verbose=false;
int juncCount=0;

std::map<std::pair<int,int>,uint64_t> donor_ends,acceptor_starts; // holds coverage data about the last/first exonic bases for the donor and acceptor sites.
std::pair<std::map<std::pair<int,int>,uint64_t>::iterator,bool> da_it;

class Intron{
public:
    Intron() = default;
    Intron(int tid,int start,int end):tid(tid),start(start),end(end){
        this->stats = std::map<int,Stats>{};
//        set_read_length(filters.max_len);
    };
    ~Intron()=default;

    void set_read_length(int read_length){
        for(auto& stat : this->stats){
            stat.second.l_lens.reserve(read_length);
            stat.second.r_lens.reserve(read_length);
            for(int i=0;i<read_length;i++){
                stat.second.l_lens[i]=0;
                stat.second.r_lens[i]=0;
            }
        }
    }

    int get_total_cov(int strand){
        return this->stats[strand].cov + this->stats[strand].cov_multi;
    }
    int get_total_cov(){
        int total_cov = 0;
        for(auto& stat : this->stats){
            total_cov += stat.second.cov + stat.second.cov_multi;
        }
        return total_cov;
    }

    int get_unique_cov(int strand){
        return this->stats[strand].cov;
    }
    int get_unique_cov(){
        int total_cov = 0;
        for(auto& stat : this->stats){
            total_cov += stat.second.cov;
        }
        return total_cov;
    }

    int get_multi_cov(int strand){
        return this->stats[strand].cov_multi;
    }
    int get_multi_cov(){
        int total_cov = 0;
        for(auto& stat : this->stats){
            total_cov += stat.second.cov_multi;
        }
        return total_cov;
    }

    double get_unique_cov_ratio(int strand){
        return (double)get_unique_cov(strand)/(double)(get_total_cov(strand));
    }
    double get_unique_cov_ratio(){
        return (double)get_unique_cov()/(double)(get_total_cov());
    }

    double get_cov_ratio(int strand){
        return (double)get_multi_cov(strand)/(double)(get_total_cov(strand));
    }
    double get_cov_ratio(){
        return (double)get_multi_cov()/(double)(get_total_cov());
    }

    double get_strand_ratio(int strand){
        return (double)get_total_cov(strand)/(double)(get_total_cov());
    }

    uint64_t get_left_cov(){
        da_it.first = donor_ends.find(std::make_pair(this->tid,this->start));
        if(da_it.first == donor_ends.end()){
            std::cerr<<"left cov does not exist...."<<std::endl;
            exit(-1);
        }
        return da_it.first->second;
    }
    uint64_t get_right_cov(){
        da_it.first = acceptor_starts.find(std::make_pair(this->tid,this->end));
        if(da_it.first == acceptor_starts.end()){
            std::cerr<<"right cov does not exist...."<<std::endl;
            exit(-1);
        }
        return da_it.first->second;
    }

    double get_left_cov_ratio(){
        return (double)get_total_cov()/(double)(get_left_cov());
    }

    double get_right_cov_ratio(){
        return (double)get_total_cov()/(double)(get_right_cov());
    }

    int get_unique_starts(int strand){
        int res = 0;
        for(int i=filters.max_len-1;i>=0;i--){
            if(this->stats[strand].l_lens[i]>0){
                res++;
            }
        }
        return res;
//        return std::set<int>(this->stats[strand].l_lens.begin(),this->stats[strand].l_lens.end() ).size();
    }
    int get_unique_starts(){
        std::set<int> res;
        for(auto& stat : stats){
            for(int i=filters.max_len-1;i>=0;i--){
                if(stat.second.l_lens[i]>0){
                    res.insert(i);
                }
            }
        }
        return res.size();
//        std::set<int> starts;
//        for(auto& stat : stats){
//            std::copy(stat.second.l_lens.begin(),stat.second.l_lens.end(),std::inserter(starts,starts.end()));
//        }
//        return starts.size();
    }

    int get_unique_ends(int strand){
        int res = 0;
        for(int i=filters.max_len-1;i>=0;i--){
            if(this->stats[strand].r_lens[i]>0){
                res++;
            }
        }
        return res;
//        return std::set<int>(this->stats[strand].r_lens.begin(),this->stats[strand].r_lens.end() ).size();
    }
    int get_unique_ends(){
        std::set<int> res;
        for(auto& stat : stats){
            for(int i=filters.max_len-1;i>=0;i--){
                if(stat.second.r_lens[i]>0){
                    res.insert(i);
                }
            }
        }
        return res.size();
//        std::set<int> ends;
//        for(auto& stat : stats){
//            std::copy(stat.second.r_lens.begin(),stat.second.r_lens.end(),std::inserter(ends,ends.end()));
//        }
//        return ends.size();
    }

    int get_max_start_dist(int strand){
        for(int i=filters.max_len-1;i>=0;i--){
            if(this->stats[strand].l_lens[i]>0){
                return i;
            }
        }
        return 0;
//        return *std::max_element(this->stats[strand].l_lens.begin(),this->stats[strand].l_lens.end());
    }
    int get_max_start_dist(){
        for(int i=filters.max_len-1;i>=0;i--){
            for(auto& stat : this->stats) {
                if(stat.second.l_lens[i]>0){
                    return i;
                }
            }
        }
        return 0;
//        std::set<int> starts;
//        for(auto& stat : stats){
//            std::copy(stat.second.l_lens.begin(),stat.second.l_lens.end(),std::inserter(starts,starts.end()));
//        }
//        return *std::max_element(starts.begin(),starts.end());
    }

    int get_max_end_dist(int strand){
        for(int i=filters.max_len-1;i>=0;i--){
            if(this->stats[strand].r_lens[i]>0){
                return i;
            }
        }
        return 0;
//        return *std::max_element(this->stats[strand].r_lens.begin(),this->stats[strand].r_lens.end());
    }
    int get_max_end_dist(){
        for(int i=filters.max_len-1;i>=0;i--){
            for(auto& stat : this->stats) {
                if(stat.second.r_lens[i]>0){
                    return i;
                }
            }
        }
        return 0;
//        std::set<int> ends;
//        for(auto& stat : stats){
//            std::copy(stat.second.r_lens.begin(),stat.second.r_lens.end(),std::inserter(ends,ends.end()));
//        }
//        return *std::max_element(ends.begin(),ends.end());
    }

    int get_expected_unique_starts(int strand){
        int max_start = get_max_start_dist(strand);
        return max_start;
    }
    int get_expected_unique_starts(){
        int max_start = get_max_start_dist();
        return max_start;
    }

    int get_expected_unique_ends(int strand){
        int max_end = get_max_end_dist(strand);
        return max_end;
    }
    int get_expected_unique_ends(){
        int max_end = get_max_end_dist();
        return max_end;
    }

    double get_frac_startsNends(int strand){ // compute fraction of expected starts and ends which are present
        int expected_starts = get_expected_unique_starts(strand);
        int expected_ends = get_expected_unique_ends(strand);
        int starts = get_unique_starts(strand);
        int ends = get_unique_ends(strand);

        return (double)(starts+ends)/(double)(expected_starts+expected_ends);
    }
    double get_frac_startsNends(){ // compute fraction of expected starts and ends which are present
        int expected_starts = get_expected_unique_starts();
        int expected_ends = get_expected_unique_ends();
        int starts = get_unique_starts();
        int ends = get_unique_ends();

        return (double)(starts+ends)/(double)(expected_starts+expected_ends);
    }

    double get_entropy_starts(int strand){
        int numlen = std::accumulate(this->stats[strand].l_lens.begin(), this->stats[strand].l_lens.end(), 0);
        double infocontent = 0;
        double expected_infocontent = 0;
        for(int i=0;i<this->stats[strand].l_lens.size();i++){
            double freq = static_cast<double>(this->stats[strand].l_lens[i])/numlen;
            if(freq==0){
                freq = 0.00000000000000001;
            }
            infocontent-= freq * log2(freq);

            double expected_freq = static_cast<double>(1)/filters.max_len;
            expected_infocontent-= expected_freq * log2(expected_freq);
        }

        return infocontent/expected_infocontent;
    }

    double get_entropy_ends(int strand){
        int numlen = std::accumulate(this->stats[strand].r_lens.begin(), this->stats[strand].r_lens.end(), 0);
        double infocontent = 0;
        double expected_infocontent = 0;
        for(int i=0;i<this->stats[strand].r_lens.size();i++){
            double freq = static_cast<double>(this->stats[strand].r_lens[i])/numlen;
            if(freq==0){
                freq = 0.00000000000000001;
            }
            infocontent-= freq * log2(freq);

            double expected_freq = static_cast<double>(1)/filters.max_len;
            expected_infocontent-= expected_freq * log2(expected_freq);
        }

        return infocontent/expected_infocontent;
    }

    double get_comb_entropy(int strand){
        double res = get_entropy_starts(strand)*get_entropy_ends(strand);
        return res;
    }

    double get_entropy_starts(){
        int numlen = 0;
        double infocontent = 0;
        double expected_infocontent = 0;

        for(auto& stat : this->stats){
            numlen += std::accumulate(stat.second.l_lens.begin(),stat.second.l_lens.end(), 0);
        }
        for(int i=0;i<filters.max_len;i++){
            double sum_strands = 0;
            for(auto& stat : this->stats){
                sum_strands += static_cast<double>(stat.second.l_lens[i]);
            }
            double freq = sum_strands/numlen;
            if(freq==0){
                freq = 0.00000000000000001;
            }
            infocontent-= freq * log2(freq);

            double expected_freq = static_cast<double>(1)/filters.max_len;
            expected_infocontent-= expected_freq * log2(expected_freq);
        }
        return infocontent/expected_infocontent;
    }

    double get_entropy_ends(){
        int numlen = 0;
        double infocontent = 0;
        double expected_infocontent = 0;

        for(auto& stat : this->stats){
            numlen += std::accumulate(stat.second.r_lens.begin(),stat.second.r_lens.end(), 0);
        }
        for(int i=0;i<filters.max_len;i++){
            double sum_strands = 0;
            for(auto& stat : this->stats){
                sum_strands += static_cast<double>(stat.second.r_lens[i]);
            }
            double freq = sum_strands/numlen;
            if(freq==0){
                freq = 0.00000000000000001;
            }
            infocontent-= freq * log2(freq);

            double expected_freq = static_cast<double>(1)/filters.max_len;
            expected_infocontent-= expected_freq * log2(expected_freq);
        }
        return infocontent/expected_infocontent;
    }

    double get_comb_entropy(){
        double res = get_entropy_starts()*get_entropy_ends();
        return res;
    }

    double get_mean_qual(int strand){
        return std::accumulate(this->stats[strand].quals.begin(),this->stats[strand].quals.end(),0.0)/(double)this->stats[strand].quals.size();
    }

    double get_mean_qual(){
        double total_sum = 0.0;
        double total_size = 0.0;
        for(auto& stat : stats){
            total_sum += std::accumulate(stat.second.quals.begin(),stat.second.quals.end(),0.0);
            total_size += stat.second.quals.size();
        }
        return total_sum/total_size;
    }

    int get_median_num_samples(int strand){
        const auto median_it = this->stats[strand].num_samples.begin() + (int)(this->stats[strand].num_samples.size() / 2);
        std::nth_element(this->stats[strand].num_samples.begin(),median_it,this->stats[strand].num_samples.end());
        return *median_it;
    }

    int get_median_num_samples(){
        std::vector<int> num_samples_all;
        for(auto& stat : stats){
            num_samples_all.insert(num_samples_all.end(),stat.second.num_samples.begin(),stat.second.num_samples.end());
        }
        const auto median_it = num_samples_all.begin() + (int)(num_samples_all.size() / 2);
        std::nth_element(num_samples_all.begin(),median_it,num_samples_all.end());
        return *median_it;
    }

    std::string get_header(){
        std::string res = "seqid\t"
                          "strand\t"
                          "start\t"
                          "end\t"
                          "med_num_samples\t"
                          "mean_qual\t"
                          "total_cov\t"
                          "unique_cov\t"
                          "unique_cov_ratio\t"
                          "strand_ratio\t"
                          "left_cov\t"
                          "right_cov\t"
                          "left_cov_ratio\t"
                          "right_cov_ratio\t"
                          "max_start\t"
                          "max_end\t"
                          "entropy_starts\t"
                          "entropy_ends\t"
                          "combined_entropy\t"
                          "unique_starts\t"
                          "unique_ends\t"
                          "frac_span";
        return res;
    }

    std::string get_tab_str_sep_strands(sam_hdr_t* al_hdr){
        std::string res = "";
        for(auto& stat : this->stats){ // for each strand
            if(stat.first==-1){
                continue;
            }
            if(stat.second.cov<=0 && stat.second.cov_multi<=0){
                continue;
            }
            res += al_hdr->target_name[this->tid];
            res += "\t" + std::string(1,(char)stat.first) +
                   "\t" + std::to_string(this->start) +
                   "\t" + std::to_string(this->end) +
                   "\t" + std::to_string(get_median_num_samples(stat.first)) +
                   "\t" + std::to_string(get_mean_qual(stat.first)) +
                   "\t" + std::to_string(get_total_cov(stat.first)) +
                   "\t" + std::to_string(get_unique_cov(stat.first)) +
                   "\t" + std::to_string(get_unique_cov_ratio(stat.first)) +
                   "\t" + std::to_string(get_strand_ratio(stat.first)) +
                   "\t" + std::to_string(get_left_cov()) +
                   "\t" + std::to_string(get_right_cov()) +
                   "\t" + std::to_string(get_left_cov_ratio()) +
                   "\t" + std::to_string(get_right_cov_ratio()) +
                   "\t" + std::to_string(get_max_start_dist(stat.first)) +
                   "\t" + std::to_string(get_max_end_dist(stat.first)) +
                   "\t" + std::to_string(get_entropy_starts(stat.first)) +
                   "\t" + std::to_string(get_entropy_ends(stat.first)) +
                   "\t" + std::to_string(get_comb_entropy(stat.first)) +
                   "\t" + std::to_string(get_unique_starts()) +
                   "\t" + std::to_string(get_unique_ends()) +
                   "\t" + std::to_string(get_frac_startsNends()) + "\n";
        }
        if(!res.empty()){
            res.pop_back();
        }
        return res;
    }

    std::string get_comb_header(){
        std::string res = "seqid\t"
                          "start\t"
                          "end\t"
                          "med_num_samples\t"
                          "mean_qual\t"
                          "total_cov\t"
                          "unique_cov\t"
                          "unique_cov_ratio\t"
                          "cov_ratio\t"
                          "left_cov\t"
                          "right_cov\t"
                          "left_cov_ratio\t"
                          "right_cov_ratio\t"
                          "max_start\t"
                          "max_end\t"
                          "entropy_starts\t"
                          "entropy_ends\t"
                          "comb_entropy\t"
                          "unique_starts\t"
                          "unique_ends\t"
                          "frac_span";
        return res;
    }

    std::string get_tab_str_combined(sam_hdr_t* al_hdr){
        std::string res = "";
        res += al_hdr->target_name[this->tid];
        res += "\t" + std::to_string(this->start) +
               "\t" + std::to_string(this->end) +
               "\t" + std::to_string(get_median_num_samples()) +
               "\t" + std::to_string(get_mean_qual()) +
               "\t" + std::to_string(get_total_cov()) +
               "\t" + std::to_string(get_unique_cov()) +
               "\t" + std::to_string(get_unique_cov_ratio()) +
               "\t" + std::to_string(get_cov_ratio()) +
               "\t" + std::to_string((int)get_left_cov()) +
               "\t" + std::to_string((int)get_right_cov()) +
               "\t" + std::to_string(get_left_cov_ratio()) +
               "\t" + std::to_string(get_right_cov_ratio()) +
               "\t" + std::to_string(get_max_start_dist()) +
               "\t" + std::to_string(get_max_end_dist()) +
               "\t" + std::to_string(get_entropy_starts()) +
               "\t" + std::to_string(get_entropy_ends()) +
               "\t" + std::to_string(get_comb_entropy()) +
               "\t" + std::to_string(get_unique_starts()) +
               "\t" + std::to_string(get_unique_ends()) +
               "\t" + std::to_string(get_frac_startsNends())+"\n";
        return res;
    }

    void add_read(bam1_t* al){ // process read to populate the stats information
        int strand = this->getXS(al);
        int multi = this->getNH(al);
        int factor = this->getYC(al);
        if(multi==1){
            this->stats[strand].cov+=(1*factor);
        }
        else{ // multimapper
            this->stats[strand].cov_multi+=(1*factor);
        }

        // look at the read and add data to the stat vectors
        this->_add_read(al);
    }
private:
    int tid=-1;
    int start=-1;
    int end=-1;

    struct Stats{
        Stats():l_lens(filters.max_len,0),r_lens(filters.max_len,0){};
        int cov=0;
        int cov_multi=0;
        int multi_weight=0; // Coverage of multimappers weighted by the number of multimappers per read
        std::vector<int> l_lens; // frequencies of lengths of the segment on the left side
        std::vector<int> r_lens; // frequencies of lengths of the segment on the right side
        std::vector<int> mismatches; // positions of mismatches relative to the intron position
        std::vector<int> quals; // quality values of reads spanning an intron
        std::vector<int> edits; // edit distances of reads
        std::vector<int> num_samples; // number of samples each read is observed in
    }; // separate statistics for two strands in case the intron is observed on both

    std::map<int,Stats> stats;

    void parse_md(std::string& md_str,std::vector<int> md_pairs){
        int pos = 0; // position on the read
    }

    void _add_read(bam1_t* al){ // populate vectors by cycling through the read positions
        int pos=al->core.pos+1; // 1-based
        int query_pos = 0; // position on the query/read

        bool first_match = true; // indicates whether the first match position has been encountered - triggers cleanup of the priority queue
        bool soft_end = false; // soft-clipped end
        int first_al_base = 0;
        int last_al_base = 0;

        int intron_coord_query = 0; // position of the intron on the query

        for (uint8_t c=0;c<al->core.n_cigar;++c){
            uint32_t *cigar_full=bam_get_cigar(al);
            int opcode=bam_cigar_op(cigar_full[c]);
            int oplen=bam_cigar_oplen(cigar_full[c]);

            switch(opcode){
                case BAM_CINS: // no change in coverage and position
                    query_pos+=oplen;
                    break;
                case BAM_CDEL: // skip to the next position - no change in coverage
                    pos+=oplen;
                    break;
                case BAM_CREF_SKIP: // skip to the next position - no change in coverage
                    // check if matches the current intron - then get left and right lengths (only aligned bases of the query - no clipping)
                    if(pos==this->start+1 && pos+oplen==this->end+1){ // current intron
                        intron_coord_query = query_pos;
                    }
                    pos+=oplen;
                    break;
                case BAM_CSOFT_CLIP: // skip to the next position - no change in coverage
                    if(!first_match){ // means this clipping is at the end of the read since a match was already found
                        soft_end = true;
                        last_al_base = query_pos;
                    }
                    query_pos+=oplen;
                    break;
                case BAM_CMATCH: // base match - increment coverage
                    if(first_match){
                        first_match=false;
                        first_al_base = query_pos;
                    }
                    for(int i=0;i<oplen;i++){ // use to look get positions of the mismatches
                        pos++;
                        query_pos++;
                    }
                    break;
                default:
                    fprintf(stderr,"ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(al));
                    exit(-1);
            }
        }
        if(!soft_end){
            last_al_base = query_pos;
        }
        int factor = getYC(al);
        int strand = getXS(al);
        int num_samples = getYX(al);
        for(int i=0;i<factor;i++){
            int r_len = last_al_base-intron_coord_query;
            int l_len = intron_coord_query-first_al_base;
            this->stats[strand].r_lens[r_len%filters.max_len]++;
            this->stats[strand].l_lens[l_len%filters.max_len]++;
            this->stats[strand].num_samples.push_back(num_samples);
            this->stats[strand].quals.push_back(al->core.qual);
        }
    }

    int getXS(bam1_t *al){
        uint8_t* ptr_xs=bam_aux_get(al,"XS");
        if(ptr_xs){
            return bam_aux2A(ptr_xs);
        }
        return -1; // -1 indicates that strand information is not available
    }

    int getNH(bam1_t *al){
        uint8_t* ptr_nh=bam_aux_get(al,"NH");
        if(ptr_nh){
            return bam_aux2i(ptr_nh);
        }
        return 1;
    }

    int getYC(bam1_t *in_rec){
        uint8_t* ptr_yc=bam_aux_get(in_rec,"YC");
        if(ptr_yc){
            return bam_aux2i(ptr_yc);
        }
        return 1;
    }

    int getYX(bam1_t *in_rec){
        uint8_t* ptr_yx=bam_aux_get(in_rec,"YX");
        if(ptr_yx){
            return bam_aux2i(ptr_yx);
        }
        return 1;
    }

    char* getMD(bam1_t *in_rec){
        uint8_t* ptr_md=bam_aux_get(in_rec,"MD");
        if(ptr_md){
            return bam_aux2Z(ptr_md);
        }
        return NULL;
    }
};

std::map<std::tuple<int,int,int>,Intron> sjs;
std::pair<std::map<std::tuple<int,int,int>,Intron>::iterator,bool> sj_it;

struct CJunc {
    int start, end;
    char strand;
    uint64_t dupcount;
    uint64_t ovh_start, ovh_end;
    CJunc(int vs=0, int ve=0, char vstrand='+', uint64_t dcount=1, uint16_t ovhs=0, uint16_t ovhe=0):
            start(vs), end(ve), strand(vstrand), dupcount(dcount),ovh_start(ovhs), ovh_end(ovhe) {}
    bool operator==(const CJunc& a) {
        return (strand==a.strand && start==a.start && end==a.end);
    }
    bool operator<(const CJunc& a) {
        if (start==a.start) return (end<a.end);
        else return (start<a.start);
    }

    void add(CJunc& j) {
        dupcount+=j.dupcount;
        if (j.ovh_start>ovh_start) ovh_start=j.ovh_start;
        if (j.ovh_end>ovh_end) ovh_end=j.ovh_end;
    }

    void write(FILE* f, const char* chr) {
        int fstart=start-ovh_start;
        int fend=end+ovh_end;
        juncCount++;
        fprintf(f, "%s\t%d\t%d\tJUNC%08d\t%u\t%c\t%d\t%d\t255,0,0,\t2\t%d,%d\t0,%d\n",
                chr, start, end, juncCount, dupcount, strand, fstart, fend,
                ovh_start, ovh_end, end-start+ovh_start);
    }
};

GArray<CJunc> junctions(128, true);

void addJunction(GSamRecord& r, int dupcount,int tid) {
    char strand = r.spliceStrand();
//    if (strand!='+' && strand!='-') return;
    for (int i=1;i<r.exons.Count();i++) {
        sj_it = sjs.insert(std::make_pair(std::make_tuple(tid,r.exons[i-1].end,r.exons[i].start-1),Intron(tid,r.exons[i-1].end,r.exons[i].start-1)));
        sj_it.first->second.add_read(r.get_b());

        // add to donor/acceptor maps
        int de = r.exons[i-1].end;
        int as = r.exons[i].start-1;
        da_it = donor_ends.insert(std::make_pair(std::make_pair(tid,r.exons[i-1].end),0));
        da_it = acceptor_starts.insert(std::make_pair(std::make_pair(tid,r.exons[i].start-1),0));

        CJunc j(r.exons[i-1].end+1, r.exons[i].start-1, strand,
                dupcount, r.exons[i-1].len(), r.exons[i].len());
        int ei;
        int r=junctions.AddIfNew(j, &ei);
        if (r==-1) {//existing junction, update
            junctions[ei].add(j);
        }
    }
}

void flushSJS(FILE* f_strand,FILE* f_comb,sam_hdr_t *al_hdr){
    std::string res_strand_sep = "";
    for(auto& sj : sjs){
        res_strand_sep = sj.second.get_tab_str_sep_strands(al_hdr);
        if(!res_strand_sep.empty()){
            fprintf(f_strand, "%s\n",res_strand_sep.c_str());
        }
        fprintf(f_comb, "%s\n",sj.second.get_tab_str_combined(al_hdr).c_str());
    }
    sjs.clear();
    donor_ends.clear();
    acceptor_starts.clear();
}

void flushJuncs(FILE* f, const char* chr) {
    for (int i=0;i<junctions.Count();i++) {
        junctions[i].write(f,chr);
    }
    junctions.Clear();
    junctions.setCapacity(128);
}

void processOptions(int argc, char* argv[]);


//b_start MUST be passed 1-based
void addCov(GSamRecord& r, int dupCount, GVec<uint64_t>& bcov, int b_start) {
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
                for(int i=0;i<oplen;i++) { // TODO: should we count deletions like this? if not counted then we end up with left/right_cov_ratio can be >1. Alternatively we could note that there is a mutation at the donor/acceptor. But this would likely require also parsing MD to also add info about SNPs
                    bcov[pos-b_start]+=dupCount;
                    pos++;
                }
                break;
            case BAM_CREF_SKIP: // skip to the next position - no change in coverage
                pos+=oplen;
                break;
            case BAM_CSOFT_CLIP:
                break;
            case BAM_CMATCH: // base match - add coverage
                for(int i=0;i<oplen;i++) {
                    bcov[pos-b_start]+=dupCount;
                    pos++;
                }
                break;
            default:
                GError("ERROR: unknown opcode: %c from read: %s",bam_cigar_opchr(opcode),bam_get_qname(in_rec));
        }
    }
}

//b_start MUST be passed 1-based
void flushCoverage(sam_hdr_t* hdr, GVec<uint64_t>& bcov,  int tid, int b_start) {
    if (tid<0 || b_start<=0) return;
    int i=0;
    b_start--; //to make it 0-based;
    while (i<bcov.Count()) {
        uint64_t icov=bcov[i];
        int j=i+1;
        while (j<bcov.Count() && icov==bcov[j]) {
            j++;
        }
        if (icov!=0)
            fprintf(coutf, "%s\t%d\t%d\t%u\n", hdr->target_name[tid], b_start+i, b_start+j, icov);

        // find those that belong to donor/acceptor sites and add coverage
        int iteration = 1;
        for(int pos=b_start+i;pos<b_start+j;pos++){
            da_it.first = donor_ends.find(std::make_pair(tid,pos+1));
            if(da_it.first != donor_ends.end()){
                da_it.first->second+=icov;
            }
            da_it.first = acceptor_starts.find(std::make_pair(tid,pos));
            if(da_it.first != acceptor_starts.end()){
                da_it.first->second+=icov;
            }
            iteration++;
        }

        i=j;
    }
}

void write_sjs(sam_hdr_t *al_hdr){
    std::fstream int_ss;
    int_ss.open(outFmt.sjs_fname,std::ios::out);

    std::fstream int_comb_ss;
    int_comb_ss.open(outFmt.sjs_combined,std::ios::out);

    int_ss<<sjs.begin()->second.get_header()<<std::endl;
    int_comb_ss<<sjs.begin()->second.get_comb_header()<<std::endl;

    std::string res_strand_sep = "";
    for(auto& sj : sjs){
        res_strand_sep = sj.second.get_tab_str_sep_strands(al_hdr);
        if(!res_strand_sep.empty()){
            int_ss<<res_strand_sep<<std::endl;
        }
        int_comb_ss<<sj.second.get_tab_str_combined(al_hdr)<<std::endl;
    }

    int_ss.close();
    int_comb_ss.close();
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
            coutf=fopen(covfname.chars(), "w");
            if (coutf==NULL) GError("Error creating file %s\n",
                                    covfname.chars());
            fprintf(coutf, "track type=bedGraph\n");
        }
    }
    if (!jfname.is_empty()) {
        joutf=fopen(jfname.chars(), "w");
        if (joutf==NULL) GError("Error creating file %s\n",
                                jfname.chars());
        fprintf(joutf, "track name=junctions\n");
    }
    if (outFmt.report_sjs) {
        sjs_strand_outf=fopen(outFmt.sjs_fname.c_str(), "w");
        if (sjs_strand_outf==NULL) GError("Error creating file %s\n",
                                          outFmt.sjs_fname.c_str());
        fprintf(sjs_strand_outf, "%s\n",Intron().get_header().c_str());

        sjs_comb_outf=fopen(outFmt.sjs_combined.c_str(), "w");
        if (sjs_comb_outf==NULL) GError("Error creating file %s\n",
                                        outFmt.sjs_combined.c_str());
        fprintf(sjs_comb_outf, "%s\n",Intron().get_comb_header().c_str());
    }
    int prev_tid=-1;
    GVec<uint64_t> bcov(4096*1024);
    int b_end=0; //bundle start, end (1-based)
    int b_start=0; //1 based
    GSamRecord brec;
    while (samreader.next(brec)) {
        int nh = brec.tag_int("NH");
        if(nh>filters.max_nh)  continue;
        if (brec.mapq()<filters.min_qual) continue;
        int endpos=brec.end;
        if (brec.refId()!=prev_tid || (int)brec.start>b_end) {
            if (coutf)
                flushCoverage(samreader.header() , bcov, prev_tid, b_start);
            if (joutf)
                flushJuncs(joutf, samreader.refName(prev_tid));
            if (outFmt.report_sjs) // +1 to account for acceptor coverage which extends 1 base past the intron
                flushSJS(sjs_strand_outf,sjs_comb_outf,samreader.header());
            b_start=brec.start;
            b_end=endpos;
            if (coutf) {
                bcov.setCount(0);
                bcov.setCount(b_end-b_start+1);
            }
            prev_tid=brec.refId();
        } else { //extending current bundle
            if (b_end<endpos) {
                b_end=endpos;
                bcov.setCount(b_end-b_start+1, (int)0);
            }
        }
        int accYC = brec.tag_int("YC", 1);
        if (coutf)
            addCov(brec, accYC, bcov, b_start);
        if (joutf && brec.exons.Count()>1) {
            addJunction(brec, accYC,prev_tid);
        }
    } //while GSamRecord emitted
    if (coutf) {
        flushCoverage(samreader.header(), bcov, prev_tid, b_start);
        if (coutf!=stdout) fclose(coutf);
    }
    if (joutf) {
        flushJuncs(joutf, samreader.refName(prev_tid));
        fclose(joutf);
    }
    if (outFmt.report_sjs) {
        flushSJS(sjs_strand_outf, sjs_comb_outf, samreader.header());
        fclose(sjs_strand_outf);
        fclose(sjs_comb_outf);
    }

//    if(outFmt.report_sjs){
//        write_sjs(samreader.header());
//    }

}// <------------------ main() end -----

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;debug;verbose;version;DVhs:c:j:b:N:Q:l:");
    args.printError(USAGE, true);
    if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
        GMessage(USAGE);
        exit(1);
    }

    if (args.getOpt('h') || args.getOpt("help")) {
        fprintf(stdout,"%s",USAGE);
        exit(0);
    }

    GStr max_nh_str=args.getOpt('N');
    if (!max_nh_str.is_empty()) {
        filters.max_nh=max_nh_str.asInt();
    }
    GStr min_qual_str=args.getOpt('Q');
    if (!min_qual_str.is_empty()) {
        filters.min_qual=min_qual_str.asInt();
    }

    GStr max_len_str=args.getOpt('l');
    if (!max_len_str.is_empty()) {
        filters.max_len=max_len_str.asInt();
    }

    debugMode=(args.getOpt("debug")!=NULL || args.getOpt('D')!=NULL);
    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);

    if (args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }
    //verbose=(args.getOpt('v')!=NULL);
    if (verbose) {
        fprintf(stderr, "Running TieCov " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    covfname=args.getOpt('c');
    jfname=args.getOpt('j');
    bfname=args.getOpt('b');
    if (args.startNonOpt()!=1) GError("Error: no alignment file given!\n");
    infname=args.nextNonOpt();

    if(args.getOpt('s')){
        outFmt.report_sjs = true;
        outFmt.sjs_fname = args.getOpt('s');
        std::string sjs_filename = args.getOpt('s');
        outFmt.sjs_combined = sjs_filename+".combined";
    }
}
