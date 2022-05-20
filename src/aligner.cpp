#include <assert.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <math.h>
#include <stdexcept>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <string>
#include <time.h>       /* time */
#include <unordered_set>
#include <utility>
#include <chrono>
#include <set>

#include "align_reads/aligner.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "edlib.h"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

#define EXTRACT nanosim_extract
#define RANDOM_SEED 422
#define OVLP_THRES 3000
#define TARGET_NUM 1000


std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace align_reads {

// ------------------------------------

struct seq_info {
    std::uint32_t id;
    bool forward;
    bool aligned;
    std::uint32_t start;
    std::uint32_t end;
    std::uint32_t left_extend;
    std::uint32_t right_extend;
    std::string contig;
    std::string name;
    std::uint32_t overlap_len;

    std::uint32_t idx_in_sequences;

    void print() {
        std::cout << "id: " << id << std::endl;
        std::cout << "forward?: " << forward << std::endl;
        std::cout << "aligned?: " << aligned << std::endl;
        std::cout << "start: " << start << std::endl;        
        std::cout << "end: " << end << std::endl;
        std::cout << "left extend: " << left_extend << std::endl;
        std::cout << "right extend: " << right_extend << std::endl;
        std::cout << "contig: " << contig << std::endl;
    }
};

struct overlap_info {
    seq_info query_info;
    seq_info target_info;

    std::uint32_t q_clip_left;
    std::uint32_t q_clip_right;
    std::uint32_t t_clip_left;
    std::uint32_t t_clip_right;
    bool is_reverse = false;
    std::uint32_t overlap_len;
    std::uint32_t score=0;

    std::uint32_t q_o_start;
    std::uint32_t q_o_end;
    std::uint32_t t_o_start;
    std::uint32_t t_o_end;

    overlap_info(seq_info q, seq_info t) : query_info(q), target_info(t) {}
};

   

struct comp {
    bool operator() (const overlap_info& lhs, const overlap_info& rhs) {
        if (lhs.query_info.id != rhs.query_info.id) return lhs.query_info.id < rhs.query_info.id;
        return lhs.target_info.id < rhs.target_info.id;
    }	
};

std::set<overlap_info, comp> all_true_overlaps;
std::set<overlap_info, comp> tps;
std::set<overlap_info, comp> fps;
std::set<overlap_info, comp> fns;

static seq_info nanosim_extract(std::string name) {
    //>NC-004354_15504699;aligned_0_R_9_1481_32
    //>NC-004354_8681245_unaligned_17149_F_0_1489_0
        
    seq_info info;
    
    std::uint32_t idx = 0;
    std::string temp;
    // extract contig name
    while (name[idx] != '_') {
        temp += name[idx++];
    }
    info.contig = temp;
    temp.clear();
    idx++;

    // extract start
    while (name[idx] != ';' && name[idx] != '_') {
        temp+= name[idx++];
    }
    info.start = std::stoi(temp);
    temp.clear();
    idx++;
    if (name[idx] == 'a') {
        info.aligned = true;
    } else {
        info.aligned = false;
    }
    while (name[idx] != '_') {
	    idx++;
    }
    idx++;
    // extract id 
    while (name[idx] != '_') {
        temp += name[idx++];
    }
    info.id = std::stoi(temp);
    temp.clear();
    idx++;
    
    // extract forward/reverse
    if (name[idx] == 'F') {
        info.forward = true;	
    } else {
	    info.forward = false;
    }
    idx += 2;

    while (name[idx] != '_') {
	    temp += name[idx++];
    }
    info.left_extend = stoi(temp);
    temp.clear();
    idx++;

    while (name[idx] != '_') {
	    temp += name[idx++];
    }
    info.end = info.start + std::stoi(temp);
    temp.clear();
    idx++;
    while (idx != name.size()) {
        temp += name[idx++];        
    }
    info.right_extend = std::stoi(temp);
    info.name = std::move(name);
    return info;    
}

static seq_info seqreq_extract(std::string name) {
    //>read=1,forward,position=19104036-19118431,length=14395,NT_033777.3
    seq_info info;

    std::uint32_t idx = 5;
    std::string temp;
    // extract read id string
    while (name[idx] != ',') {
        temp += name[idx++];
    }
    info.id = std::stoi(temp);
    temp.clear();
    idx++;

    // extract forward or reverse
    while (name[idx] != ',') {
        temp+= name[idx++];
    }
    if (temp.compare("forward") == 0) {
	info.forward = true;	
    } else {
	info.forward = false;
    }
    temp.clear();
    idx += 10;

    // extract start 
    while (name[idx] != '-') {
        temp += name[idx++];
    }
    info.start = std::stoi(temp);
    temp.clear();
    idx++;

    // extract end
    while (name[idx] != ',') {
        temp += name[idx++];
    }
    info.end = std::stoi(temp);
    temp.clear();
    idx++;
    while (name[idx] != ',') {
	idx++;
    }
    idx++;

    // extract contig
    while (idx < name.size()) {
	temp += name[idx++];
    }
    info.contig = std::move(temp);

    return info;    


}

std::uint32_t overlap_len_sub(seq_info query_info, seq_info target_info) {
	if (query_info.contig != target_info.contig) return 0;     
    if (query_info.start >= target_info.start && query_info.start < target_info.end) {
        return query_info.end < target_info.end ? query_info.end - query_info.start : target_info.end - query_info.start;
	       
	} else if (query_info.end <= target_info.end && query_info.end > target_info.start) {
	    return query_info.start < target_info.start ? query_info.end - target_info.start : query_info.end - query_info.start;
	}
	return 0;
}

std::uint32_t overlap_len(seq_info query_info, seq_info target_info) {
    return std::max(overlap_len_sub(query_info, target_info), overlap_len_sub(target_info, query_info));	
}




//------------------------------------------------

void Aligner::run() {
    find_true_overlaps();   
    within_each();
    //find_RAM_overlaps(); 
    //true_positives();
    //false_positives();
    //false_negatives();
}


void Aligner::find_true_overlaps() {
    for (auto i = 0; i < TARGET_NUM; i++) {        

        auto target_info = EXTRACT(sequences[i]->name);
        for (auto j = 0; j < sequences.size(); j++) {
            if (j == i) continue;
            auto query_info = EXTRACT(sequences[j]->name);                       
            std::uint32_t len = overlap_len(query_info, target_info);

            if (len > OVLP_THRES) {
                query_info.idx_in_sequences = j;
                target_info.idx_in_sequences = i;
                overlap_info o_info {query_info, target_info};
                o_info.overlap_len = len;                
                all_true_overlaps.insert(o_info);
            }   
        } 
    }
}


void Aligner::within_each() {
    

    for (auto i = 0; i < TARGET_NUM; i++) {
        auto& target = sequences[i];
        std::uint32_t target_id = target->id;
        auto target_info = EXTRACT(target->name);
        target_info.idx_in_sequences = id_to_pos_index[target_id];
        std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, true);

        std::uint32_t num_tp = 0;
        std::uint32_t num_fp = 0;
        
        double total_norm_dist_tp = 0;
        double total_norm_score_whole_tp = 0;
        double total_norm_score_part_tp = 0;


        double total_norm_dist_fp = 0;
        double total_norm_score_whole_fp = 0;
        double total_norm_score_part_fp = 0;

        for (auto& o : overlaps) {
            auto& s = sequences[id_to_pos_index[o.rhs_id]];
            auto query_info = EXTRACT(s->name);
            query_info.idx_in_sequences = id_to_pos_index[o.rhs_id];            
            overlap_info o_info {query_info, target_info};
 

            std::uint32_t t_begin = o.lhs_begin;
            std::uint32_t t_end = o.lhs_end;
            std::uint32_t q_begin = o.rhs_begin;
            std::uint32_t q_end = o.rhs_end;

            if (!o.strand) {
                s->ReverseAndComplement();
                o_info.is_reverse = true;
                q_end = s->inflated_len - o.rhs_begin;
                q_begin = s->inflated_len - o.rhs_end;
            }
            int protrude_left = q_begin - t_begin;
            int protrude_right = (s->inflated_len - q_end) - (target->inflated_len - t_end);

            std::uint32_t q_clip_left;
            std::uint32_t q_clip_right;
            std::uint32_t t_clip_left;
            std::uint32_t t_clip_right;
            if (protrude_left > 0) {
                q_clip_left = protrude_left;
                t_clip_left = 0;
            } else {
                q_clip_left = 0;
                t_clip_left = -protrude_left;
            }
            if (protrude_right > 0) {
                q_clip_right = protrude_right;
                t_clip_right = 0;
            } else {
                q_clip_right = 0;
                t_clip_right = -protrude_right;
            }
  
                       
            o_info.q_clip_left = q_clip_left;
            o_info.q_clip_right = q_clip_right;
            o_info.t_clip_left = t_clip_left;
            o_info.t_clip_right = t_clip_right;

            o_info.q_o_start = o.rhs_begin;
            o_info.q_o_end = o.rhs_end;
            o_info.t_o_start = o.lhs_begin;
            o_info.t_o_end = o.lhs_end;
            o_info.score = o.score;

            std::uint32_t query_string_len = s->inflated_len - q_clip_left - q_clip_right;
            std::uint32_t target_string_len = target->inflated_len - t_clip_left - t_clip_right;
            std::string query_string = s->InflateData(q_clip_left, query_string_len);
            std::string target_string = target->InflateData(t_clip_left, target_string_len);
            EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                target_string.c_str(), target_string.size(),
                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            if (!o.strand) {
                s->ReverseAndComplement();
            }

            if (overlap_len(query_info, target_info) > OVLP_THRES) {
                tps.insert(o_info);
                total_norm_dist_tp += (double) result.editDistance/query_string_len;
                total_norm_score_whole_tp += (double) o.score/std::max(query_string_len, target_string_len);
                total_norm_score_part_tp += (double) o.score/std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.rhs_begin);
                num_tp++;


            } else {
                fps.insert(o_info);
                total_norm_dist_fp += (double) result.editDistance/query_string_len;
                total_norm_score_whole_fp += (double) o.score/std::max(query_string_len, target_string_len);
                total_norm_score_part_fp += (double) o.score/std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.rhs_begin);
                num_fp++;

            }


            edlibFreeAlignResult(result);
        }
        std::cout << "-------------------------------------------" << std::endl;
        std::cout << "num tp: " << num_tp << std::endl;    
        std::cout << "avg norm dist tp : " << total_norm_dist_tp/num_tp << std::endl;
        std::cout << "avg norm score whole tp: " << total_norm_score_whole_tp/num_tp << std::endl;
        std::cout << "avg norm score part tp: " << total_norm_score_part_tp/num_tp << std::endl;
   
        std::cout << "-----------------------" << std::endl;
        std::cout << "num fp: " << num_fp << std::endl;    
        std::cout << "avg norm dist fp : " << total_norm_dist_fp/num_fp << std::endl;
        std::cout << "avg norm score whole fp: " << total_norm_score_whole_fp/num_fp << std::endl;
        std::cout << "avg norm score part fp: " << total_norm_score_part_fp/num_fp << std::endl;
        std::cout << "-------------------------------------------" << std::endl;

   }

    for (auto& o_info: all_true_overlaps) {
        if (tps.find(o_info) == tps.end()) {
            fns.insert(o_info);
        }
    }

    std::cout << "precision: " << ((double) tps.size() / (tps.size() + fps.size())) << std::endl;
    std::cout << "recall: " << ((double) tps.size() / all_true_overlaps.size()) << std::endl;
}




void Aligner::find_RAM_overlaps() {

    

    for (auto i = 0; i < TARGET_NUM; i++) {
        auto& target = sequences[i];
        std::uint32_t target_id = target->id;
        auto target_info = EXTRACT(target->name);
        target_info.idx_in_sequences = id_to_pos_index[target_id];
        std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, true);          

        for (auto& o : overlaps) {
            auto& s = sequences[id_to_pos_index[o.rhs_id]];
            auto query_info = EXTRACT(s->name);
            query_info.idx_in_sequences = id_to_pos_index[o.rhs_id];            
            overlap_info o_info {query_info, target_info};
 

            std::uint32_t t_begin = o.lhs_begin;
            std::uint32_t t_end = o.lhs_end;
            std::uint32_t q_begin = o.rhs_begin;
            std::uint32_t q_end = o.rhs_end;

            if (!o.strand) {
                o_info.is_reverse = true;
                q_end = s->inflated_len - o.rhs_begin;
                q_begin = s->inflated_len - o.rhs_end;
            }
            int protrude_left = q_begin - t_begin;
            int protrude_right = (s->inflated_len - q_end) - (target->inflated_len - t_end);

            std::uint32_t q_clip_left;
            std::uint32_t q_clip_right;
            std::uint32_t t_clip_left;
            std::uint32_t t_clip_right;
            if (protrude_left > 0) {
                q_clip_left = protrude_left;
                t_clip_left = 0;
            } else {
                q_clip_left = 0;
                t_clip_left = -protrude_left;
            }
            if (protrude_right > 0) {
                q_clip_right = protrude_right;
                t_clip_right = 0;
            } else {
                q_clip_right = 0;
                t_clip_right = -protrude_right;
            }
  
                       
            o_info.q_clip_left = q_clip_left;
            o_info.q_clip_right = q_clip_right;
            o_info.t_clip_left = t_clip_left;
            o_info.t_clip_right = t_clip_right;

            o_info.q_o_start = o.rhs_begin;
            o_info.q_o_end = o.rhs_end;
            o_info.t_o_start = o.lhs_begin;
            o_info.t_o_end = o.lhs_end;
            o_info.score = o.score;

            if (overlap_len(query_info, target_info) > OVLP_THRES) {
                tps.insert(o_info);
            } else {
                fps.insert(o_info);
            }
        }
    }

    for (auto& o_info: all_true_overlaps) {
        if (tps.find(o_info) == tps.end()) {
            fns.insert(o_info);
        }
    }

    std::cout << "precision: " << ((double) tps.size() / (tps.size() + fps.size())) << std::endl;
    std::cout << "recall: " << ((double) tps.size() / all_true_overlaps.size()) << std::endl;
}



void Aligner::false_negatives() {
    std::cout << "----false negatives----" << std::endl; 
    // true positives
    for (auto& o_info : fns) {
        if (!o_info.query_info.aligned || !o_info.target_info.aligned) {
            continue;
        }
         
        auto& query_sequence = sequences[o_info.query_info.idx_in_sequences];
        auto& target_sequence = sequences[o_info.target_info.idx_in_sequences];


        if (!o_info.query_info.forward) query_sequence->ReverseAndComplement();
        if (!o_info.target_info.forward) target_sequence->ReverseAndComplement();
        std::uint32_t q_start = o_info.query_info.start;
        std::uint32_t q_end = o_info.query_info.end;
        std::uint32_t t_start = o_info.target_info.start;
        std::uint32_t t_end = o_info.target_info.end;

        int q_clip_left = t_start - q_start;
        int q_clip_right = q_end - t_end;
        int t_clip_left;
        int t_clip_right;
        if (q_clip_left >= 0) {
            t_clip_left = 0;
        } else {
            t_clip_left = - q_clip_left;
            q_clip_left = 0;
        }
        if (q_clip_right >= 0) {
            t_clip_right = 0;
        } else {
            t_clip_right = -q_clip_right;
            q_clip_right = 0;
        }
       
        q_clip_left += o_info.query_info.left_extend;
        q_clip_right += o_info.query_info.right_extend;
        t_clip_left += o_info.target_info.left_extend;
        t_clip_right += o_info.target_info.right_extend;
       
        
        std::uint32_t ovlp_len = overlap_len(o_info.query_info, o_info.target_info);
        int query_string_len = query_sequence->inflated_len - q_clip_left - q_clip_right;
        int target_string_len = target_sequence->inflated_len - t_clip_left - t_clip_right;
        if (query_string_len == 0 || target_string_len==0 || ovlp_len < 0.1 * std::max(query_sequence->inflated_len, target_sequence->inflated_len)) {
            continue;
        }
  
        std::cout << "query name: " << o_info.query_info.name << std::endl;
        std::cout << "target name: " << o_info.target_info.name << std::endl;

        std::cout << "pre-clip q-len: " << query_sequence->inflated_len << " pre-clip t-len: " << target_sequence->inflated_len << std::endl;
        std::cout << "q_start: " << o_info.query_info.start << " query_end: " << o_info.query_info.end << std::endl;
        std::cout << "query_clip_left: " << q_clip_left << " query_clip_right: " << q_clip_right << std::endl;
        std::cout << "t_start: " << o_info.target_info.start << " target_end: " << o_info.target_info.end << std::endl;
        std::cout << "t_clip_left: " << t_clip_left << " target_clip_right: " << t_clip_right << std::endl;
        std::cout << "overlap len: " << ovlp_len << std::endl;
        
      

        std::string query_string = query_sequence->InflateData(q_clip_left, query_string_len);
        std::string target_string = target_sequence->InflateData(t_clip_left, target_string_len);
       


        if (!o_info.query_info.forward) query_sequence->ReverseAndComplement();
        if (!o_info.target_info.forward) target_sequence->ReverseAndComplement();

        std::vector<std::string> overlapping_reads;
       

        EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                        target_string.c_str(), target_string.size(),
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

        overlapping_reads.push_back(std::move(query_string));
        std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
        clips.emplace_back(0, 0);

        std::vector<EdlibAlignResult> rs = {result};

        auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
        
        alignment.print();

                     
    }
}




void Aligner::true_positives() {

    std::cout << "----true positives----" << std::endl; 
    
    double total_norm_dist = 0;
    double total_norm_score_whole = 0;
    double total_norm_score_part = 0;

    // true positives
    for (auto& o_info : tps) {
         
         auto& query_sequence = sequences[o_info.query_info.idx_in_sequences];
         auto& target_sequence = sequences[o_info.target_info.idx_in_sequences];
         
         if (o_info.is_reverse) query_sequence->ReverseAndComplement();
         std::uint32_t query_string_len = query_sequence->inflated_len - o_info.q_clip_left - o_info.q_clip_right;
         std::uint32_t target_string_len = target_sequence->inflated_len - o_info.t_clip_left - o_info.t_clip_right;
         std::string target_string = target_sequence-> InflateData(o_info.t_clip_left, target_string_len);
         std::string query_string = query_sequence->InflateData(o_info.q_clip_left, query_string_len);

         EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
             target_string.c_str(), target_string.size(),
             edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));     

             total_norm_dist += (double) result.editDistance/query_string_len;         
             total_norm_score_whole = (double) o_info.score / std::max(query_string_len, target_string_len);
             total_norm_score_part = (double) o_info.score / std::max(o_info.t_o_end - o_info.t_o_start, o_info.q_o_end - o_info.q_o_start); 

             std::vector<std::string> overlapping_reads;
             overlapping_reads.push_back(std::move(query_string));
             std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
             clips.emplace_back(0, 0);
             std::vector<EdlibAlignResult> rs = {result};
             auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
             std::cout << "query name: " << o_info.query_info.name << std::endl;
             std::cout << "target name: " << o_info.target_info.name << std::endl;
             std::cout << "query_start: " << o_info.query_info.start << " query_end: " << o_info.query_info.end << std::endl;
             std::cout << "q_clip_left: " << o_info.q_clip_left << " query_clip_right: " << o_info.q_clip_right << std::endl;
             std::cout << "t_start: " << o_info.target_info.start << " target_end: " << o_info.target_info.end << std::endl;
             std::cout << "t_clip_left: " << o_info.t_clip_left << " target_clip_right: " << o_info.t_clip_right << std::endl;
             std::uint32_t ovlp_len = overlap_len(o_info.query_info, o_info.target_info);
             std::cout << "overlap len: " << ovlp_len << std::endl;
             //std::cout << "q_o_left: " << o_info.q_o_left << " query_o_right: " << o_info.q_o_right << std::endl;
             //std::cout << "t_o_left: " << o_info.t_o_left << " target_o_right: " << o_info.t_o_right << std::endl;
             alignment.print();
             if (o_info.is_reverse) query_sequence->ReverseAndComplement();
    }

    std::cout << "avg norm dist: " << total_norm_dist / fps.size() << std::endl;
    std::cout << "avg norm score whole: " << total_norm_score_whole / fps.size() << std::endl;
    std::cout << "avg norm score part: " << total_norm_score_part / fps.size() << std::endl;
    std::cout << "----true positives----" << std::endl; 

}

void Aligner::false_positives() {

    std::cout << "----false positives----" << std::endl; 

    double total_norm_dist = 0;
    double total_norm_score_whole = 0;
    double total_norm_score_part = 0;

    // true positives
    for (auto& o_info : fps) {
         
         auto& query_sequence = sequences[o_info.query_info.idx_in_sequences];
         auto& target_sequence = sequences[o_info.target_info.idx_in_sequences];
         
         if (o_info.is_reverse) query_sequence->ReverseAndComplement();
         
         std::uint32_t query_string_len = query_sequence->inflated_len - o_info.q_clip_left - o_info.q_clip_right;
         std::uint32_t target_string_len = target_sequence->inflated_len - o_info.t_clip_left - o_info.t_clip_right;
         std::string target_string = target_sequence-> InflateData(o_info.t_clip_left, target_string_len);
         std::string query_string = query_sequence->InflateData(o_info.q_clip_left, query_string_len);
         EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
             target_string.c_str(), target_string.size(),
             edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

             total_norm_dist += (double) result.editDistance/query_string_len;         
             total_norm_score_whole = (double) o_info.score / std::max(query_string_len, target_string_len);
             total_norm_score_part = (double) o_info.score / std::max(o_info.t_o_end - o_info.t_o_start, o_info.q_o_end - o_info.q_o_start); 


             std::vector<std::string> overlapping_reads;
             overlapping_reads.push_back(std::move(query_string));
             std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
             clips.emplace_back(0, 0);
             std::vector<EdlibAlignResult> rs = {result};
             auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
             std::cout << "query name: " << o_info.query_info.name << std::endl;
             std::cout << "target name: " << o_info.target_info.name << std::endl;
             std::cout << "query_start: " << o_info.query_info.start << " query_end: " << o_info.query_info.end << std::endl;
             std::cout << "q_clip_left: " << o_info.q_clip_left << " query_clip_right: " << o_info.q_clip_right << std::endl;
             std::cout << "t_start: " << o_info.target_info.start << " target_end: " << o_info.target_info.end << std::endl;
             std::cout << "t_clip_left: " << o_info.t_clip_left << " target_clip_right: " << o_info.t_clip_right << std::endl;
             std::uint32_t ovlp_len = overlap_len(o_info.query_info, o_info.target_info);
             std::cout << "overlap len: " << ovlp_len << std::endl;


             //std::cout << "q_o_left: " << o_info.q_o_left << " query_o_right: " << o_info.q_o_right << std::endl;
             //std::cout << "t_o_left: " << o_info.t_o_left << " target_o_right: " << o_info.t_o_right << std::endl;
             alignment.print();
             if (o_info.is_reverse) query_sequence->ReverseAndComplement();
    }

    std::cout << "avg norm dist: " << total_norm_dist / fps.size() << std::endl;
    std::cout << "avg norm score whole: " << total_norm_score_whole / fps.size() << std::endl;
    std::cout << "avg norm score part: " << total_norm_score_part / fps.size() << std::endl;
    std::cout << "----false positives----" << std::endl; 

}



Aligner::Aligner(const char** sequences_paths, std::shared_ptr<thread_pool::ThreadPool>& pool, std::uint8_t kmer_len, 
    std::uint8_t window_len, double freq) 
    : pool(pool), minimizer_engine(pool, kmer_len, window_len, 500, 4, 100, 500) {
    
	
    //srand (time(NULL));
    srand(RANDOM_SEED);
    auto is_suffix = [] (const std::string& s, const std::string& suff) {
        return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
    };    
    for (std::uint8_t i = 0; i < 2; i++) {
        const char* path = sequences_paths[i];

        if (path == nullptr) break;
        std::string seq_path {path};
        if (is_suffix(seq_path, ".fasta")) {
            try { 
                auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(seq_path);
                auto s = p->Parse(-1);
                sequences.insert(sequences.end(), std::make_move_iterator(s.begin()), std::make_move_iterator(s.end()));
            } catch (const std::invalid_argument& exception) {
                std::cerr << exception.what() << std::endl;
            }            
                  
        } else if (is_suffix(seq_path, ".fastq")) {
            try { 
                auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(seq_path);
                auto s = p->Parse(-1);
                sequences.insert(sequences.end(), std::make_move_iterator(s.begin()), std::make_move_iterator(s.end()));
            } catch (const std::invalid_argument& exception) {
                std::cerr << exception.what() << std::endl;
            }            

        } else {
            throw std::invalid_argument("[align_reads::Aligner::Aligner] Error: Invalid sequences file format.");
        }
    }
    auto long_seq_first = [] (const std::unique_ptr<biosoup::NucleicAcid>& s1, const std::unique_ptr<biosoup::NucleicAcid>& s2) {
        return s1->inflated_len > s2->inflated_len;        
    };
    

    
    id_to_pos_index.resize(sequences.size());
    
    std::random_shuffle(sequences.begin(), sequences.end());
    
 
    std::uint32_t pos_index = 0;

    for (auto& s: sequences) {
        id_to_pos_index[s->id] = pos_index++; 
    }
	
       
    minimizer_engine.Minimize(sequences.begin(), sequences.end(), true);
    minimizer_engine.Filter(freq);
     
	
}

Aligner::align_result Aligner::align_to_target_clip(std::vector<std::string>& queries, std::string& target,
    std::vector<std::pair<std::uint32_t, std::uint32_t>>& clips, std::vector<std::vector<std::uint32_t>>& inserters,
    std::vector<std::uint32_t>& ins_at_least2, bool has_hap, std::vector<EdlibAlignResult>& edlib_results) {  
    align_result result;
    // Fill in the bases from the target into the first row
    result.target_columns.resize(target.size());        
    for (std::uint32_t i = 0; i < result.target_columns.size(); i++) { // for each column
        auto& column = result.target_columns[i];
        column.resize(queries.size() + 1, '_');
        column[0] = target[i]; // the first item of each row - makes up the first row
    }
    // Allowance for insertion columns after each target position
    result.ins_columns.resize(target.size());
    
    // Initialize width of the alignment to the length of the target 
    result.width = target.size(); // will be incremented along with adding insertion columns
    
    // Initialize the record of who is inserting at each position
    
    
    // Aligning queries to target and filling up/adding columns
    for (std::uint32_t k = 0; k < queries.size(); k++) {   
        auto& clip = clips[k];
        std::uint32_t left_clip = clip.first;
        std::uint32_t right_clip = clip.second;        
        auto& query = queries[k];
        EdlibAlignResult& align = edlib_results[k];

		
		/*edlibAlign(query.c_str(), query.size(),
            target.c_str() + left_clip, target.size() - left_clip - right_clip,
            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));*/
            
        std::uint32_t next_query_index = 0; // upcoming query position
        std::uint32_t next_target_index = align.startLocations[0] + left_clip; // upcoming target position
        std::uint32_t next_ins_index = 0; // upcoming ins index for the current segment
        std::uint32_t align_start = next_target_index;
        std::uint32_t align_end = static_cast<std::uint32_t>(align.endLocations[0] + left_clip);
       
        for (int i = 0; i < align.alignmentLength; i++) {
            
            switch (align.alignment[i]) {
                case 0: { // match

                    //place matching query base from the query into column
                   
                    result.target_columns[next_target_index++][k+1] = query[next_query_index++]; 

                    next_ins_index = 0;
                    break;
                }
                case 1: { // ins
                    std::uint32_t ins_columns_index = next_target_index - 1;
                   
                    // insertions in the query to the sides of target is ignored
                    if (next_target_index == 0) {
                        next_query_index++;
                        continue;
                    } else if (next_target_index == target.size()) {
                        i = align.alignmentLength;  // if at ins of the query to the right of target, terminate loop early
                        continue;                       
                    }    
                    
                                        
                    // check if we have enough columns ready for this target position
                    auto& ins_columns = result.ins_columns[ins_columns_index];
                    if (ins_columns.size() < next_ins_index + 1) {
                        // ignore extra ins column due to hap sequence
                        if (has_hap && k >= queries.size() - 2) {
                            next_query_index++;
                            continue;
                        }
                        // if not enough, create new column
                        ins_columns.resize(next_ins_index + 1);
                        ins_columns[next_ins_index].resize(queries.size() + 1, '_');
                        result.width++;
                    }
                    
                    // Record existence of ins at certain positions
                    if (next_ins_index == 0) {
                        auto& v = inserters[ins_columns_index]; 
                        v.push_back(k);
                        if (v.size() == 2) {
                            ins_at_least2.push_back(ins_columns_index);
                        }
                    }
                    
                    ins_columns[next_ins_index++][k+1] = query[next_query_index++];
                    break;
                }
                case 2: { // del
                    next_target_index++;
                    next_ins_index = 0;
                    break;
                }
                case 3: { // mismatch
                    // place mismatching query base into column
                    result.target_columns[next_target_index++][k+1] = query[next_query_index++]; 
                    next_ins_index = 0;
                    break;
                }
                default: {
                    std::cerr << "Unknown alignment result by edlib!" << std::endl;
                }        
            }            
        }
        
        edlibFreeAlignResult(align); 
    }
    
    return result;
}

Aligner::align_result Aligner::align_to_target_no_clip(std::vector<std::string>& queries, std::string& target, bool has_hap) {
     
    align_result result;
    // Fill in the bases from the target into the first row
    result.target_columns.resize(target.size());        
    for (std::uint32_t i = 0; i < result.target_columns.size(); i++) { // for each column
        auto& column = result.target_columns[i];
        column.resize(queries.size() + 1, '_');
        column[0] = target[i]; // the first item of each row - makes up the first row
    }

    // Allowance for insertion columns after each target position
    result.ins_columns.reserve(target.size() + 1);
    result.ins_columns.resize(target.size());
    
    // Initialize width of the alignment to the length of the target 
    result.width = target.size(); // will be incremented along with adding insertion columns
    
    // Aligning queries to target and filling up/adding columns
    for (std::uint32_t k = 0; k < queries.size(); k++) {      

        auto& query = queries[k];
        EdlibAlignResult align = edlibAlign(query.c_str(), query.size(),
            target.c_str(), target.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
        std::uint32_t next_query_index = 0; // upcoming query position
        std::uint32_t next_target_index = align.startLocations[0]; // upcoming target position
        std::uint32_t next_ins_index = 0; // upcoming ins index for the current segment
        std::uint32_t align_start = next_target_index;
        std::uint32_t align_end = static_cast<std::uint32_t>(align.endLocations[0]);
        bool contained = true;
        
        for (int i = 0; i < align.alignmentLength; i++) {
            
         
            switch (align.alignment[i]) {
                case 0: { // match
 
                    //place matching query base from the query into column
                    
                    result.target_columns[next_target_index++][k+1] = query[next_query_index++]; 

                    next_ins_index = 0;
                    break;
                }
                case 1: { // ins
                   
                    std::uint32_t ins_columns_index = next_target_index - 1;
                    
                    //Ins to left of target. 
                    if (next_target_index == 0) {
                        contained = false;
                        if (has_hap && k >= queries.size() - 2) {
                            next_query_index++;
                            continue;
                        }
                        result.ins_columns.resize(target.size() + 1);
                        ins_columns_index = target.size();        
                    }
                                           
                    // check if we have enough columns ready for this target position
                    auto& ins_columns = result.ins_columns[ins_columns_index];
                    if (ins_columns.size() < next_ins_index + 1) {
                        if (has_hap && k >= queries.size() - 2) {
                            next_query_index++;
                            continue;
                        }
                        // if not enough, create new column
                        ins_columns.resize(next_ins_index + 1);
                        ins_columns[next_ins_index].resize(queries.size() + 1, '_');
                        result.width++;
                    } 
                    ins_columns[next_ins_index++][k+1] = query[next_query_index++];
                    break;
                }
                case 2: { // del
                    
                    next_target_index++;
                    next_ins_index = 0;
                    break;
                }
                case 3: { // mismatch
                   
                    // place mismatching query base into column
                    result.target_columns[next_target_index++][k+1] = query[next_query_index++]; 
                    next_ins_index = 0;
                    break;
                }
                default: {
                    std::cerr << "Unknown alignment result by edlib!" << std::endl;
                }        
            }            
        }
        
        
        edlibFreeAlignResult(align); 
    }
    
    return result;
}




Aligner::align_result Aligner::multi_align(std::vector<std::string>& queries, std::string& target, 
    std::vector<std::pair<std::uint32_t, std::uint32_t>>& clips, std::vector<EdlibAlignResult>& edlib_results, bool has_hap) {
    std::vector<std::vector<std::uint32_t>> inserters; //for each position, contains positional index of queries that has ins there    
    inserters.resize(target.size());
    std::vector<std::uint32_t> ins_at_least2;    // target positions where there are at least 2 reads with ins
    
    align_result result = align_to_target_clip(queries, target, clips, inserters, ins_at_least2, has_hap, edlib_results);
    
    for (auto& ins_pos : ins_at_least2) { // for all positions that have at least two queries with insertions, 
        // extract the ins segments, and record the longest one  
        std::uint32_t position_index_longest = 0; // index in ins_segments
        std::uint32_t longest_ins_len = 0;
        
        std::vector<std::string> ins_segments;
        auto inserters_at_pos = inserters[ins_pos]; // all the inserters here
        ins_segments.resize(inserters_at_pos.size());
        auto& ins_columns = result.ins_columns[ins_pos];
        for (auto& column : ins_columns) { // for each column, extract the non gap bases for the queries with ins there
            for (std::uint32_t i = 0 ; i < inserters_at_pos.size(); i++) {
                char& base = column[inserters_at_pos[i]+1];
                auto& segment = ins_segments[i];
                if (base != '_') {
                    segment.push_back(base);
                    if (segment.size() > longest_ins_len && (!has_hap || inserters_at_pos[i] < queries.size() - 2)) {
                        position_index_longest = i;
                        longest_ins_len = segment.size();
                    }
                }
            }            
        }
        
        // align the rest to the longest segment
        std::uint32_t query_position_longest_ins = inserters_at_pos[position_index_longest];
        inserters_at_pos.erase(inserters_at_pos.begin() + position_index_longest); //exclude longest from list of queries 
        std::string longest_ins_segment = std::move(ins_segments[position_index_longest]);
        ins_segments.erase(ins_segments.begin() + position_index_longest);
        auto sub_result = align_to_target_no_clip(ins_segments, longest_ins_segment, has_hap);
        
        std::uint32_t t = ins_columns.size(); 

        // get more ins columns for a target position if needed
        ins_columns.resize(sub_result.width);
        result.width += sub_result.width - t;
        for (; t < ins_columns.size(); t++) {
            ins_columns[t].resize(queries.size() + 1 ,'_');	
        }

        std::uint32_t to_column_id = 0;
        
        // if there is ins to the left of target...
        if (sub_result.ins_columns.size() > sub_result.target_columns.size()) {
            auto& sub_ins_columns = sub_result.ins_columns[sub_result.target_columns.size()];
            for (std::uint32_t j = 0; j < sub_ins_columns.size(); j++) {
                auto& to_column = ins_columns[to_column_id++];
                auto& from_column = sub_ins_columns[j];
                to_column[query_position_longest_ins + 1] = from_column[0];    
                for (std::uint32_t j = 0; j < inserters_at_pos.size(); j++) {
                    to_column[inserters_at_pos[j] + 1] = from_column[j + 1];               
                }    
            }
            
        }
        
        for(std::uint32_t i = 0; i < sub_result.target_columns.size(); i++) {
            auto& to_column = ins_columns[to_column_id++];
            auto& from_column = sub_result.target_columns[i];
            
            // traversing the from columns 
            to_column[query_position_longest_ins + 1] = from_column[0]; // the target is the top row of the alignment
            
            for (std::uint32_t j = 0; j < inserters_at_pos.size(); j++) {
                    to_column[inserters_at_pos[j] + 1] = from_column[j + 1];             
                }

                auto& sub_ins_columns = sub_result.ins_columns[i];
                
            for (std::uint32_t j = 0; j < sub_ins_columns.size(); j++) {
                    auto& to_column = ins_columns[to_column_id++];
                    
            auto& from_column = sub_ins_columns[j];
                    
            to_column[query_position_longest_ins + 1] = from_column[0];    
                    
            for (std::uint32_t j = 0; j < inserters_at_pos.size(); j++) {
                        to_column[inserters_at_pos[j] + 1] = from_column[j + 1];               
                    }    
               
            }
            	    
        }
		
    }

    
    return result;
}



void Aligner::align_result::print() {
    std::uint32_t block_width = 100;
    std::vector<std::vector<std::string>> output_blocks;
    std::uint32_t num_rows = this->target_columns[0].size() * 2 -1;
    std::uint32_t num_blocks = this->width / block_width;
    num_blocks = this->width % block_width == 0 ? num_blocks : num_blocks + 1;
    output_blocks.resize(num_blocks);
    for (auto& block : output_blocks) {
        block.resize(num_rows);
        for (auto& row : block) {
            row.reserve(width);
        }
    }
    std::vector<std::uint32_t> start_of_blocks;
    start_of_blocks.reserve(num_blocks);
    start_of_blocks.push_back(0);
    std::uint32_t col_index = 0;
    std::uint32_t current_block_index = 0;
    if (this->ins_columns.size() > this->target_columns.size()) {
        
        auto& ins_cols = this->ins_columns[this->target_columns.size()];
        for (std::uint32_t j = 0; j < ins_cols.size() ; j++) { // for each ins column here
            auto& ins_column = ins_cols[j];
	    std::uint32_t block_index = col_index++/block_width;
            auto& block = output_blocks[block_index];
	    if (block_index > current_block_index) {
	        current_block_index++;
		start_of_blocks.push_back(0);	
	    }
            for (std::uint32_t k = 0; k < num_rows; k++) { // move down the column
                if (k%2==0) { 
                    block[k].push_back(ins_column[k/2]);
                    if (k != 0) {
                        if (block[k].back() != '_' && block[k-2].back() != '_' && block[k].back() != block[k-2].back()) {
                            block[k-1].push_back('!');
                        } else {
                            block[k-1].push_back(' ');
                        }
                    }
                }
            }
        }     
    }
    for (std::uint32_t i = 0; i < this->target_columns.size(); i++) { // for each target position        
	auto& column = this->target_columns[i];
	
	std::uint32_t block_index = col_index++/block_width;
        auto& block = output_blocks[block_index];
        if (block_index > current_block_index) {
	        current_block_index++;
		start_of_blocks.push_back(i);	
	}

	for (std::uint32_t k = 0; k < num_rows; k++) { // move down the column
	    
            if (k%2==0) {
                block[k].push_back(column[k/2]);

                if (k != 0) {
			
                    if (block[k].back() != '_' && block[k-2].back() != '_' && block[k].back() != block[k-2].back()) {
                        
			block[k-1].push_back('!');
                    } else {
			    
                        block[k-1].push_back(' ');
                    }

                }
            }        
        }

        auto& ins_cols = this->ins_columns[i];
        
	for (std::uint32_t j = 0; j < ins_cols.size(); j++) { // for each ins column here
	    std::uint32_t block_index = col_index++/block_width;	
            auto& block = output_blocks[block_index];
            if (block_index > current_block_index) {
	        current_block_index++;
                start_of_blocks.push_back(i);	
	    }

            auto& ins_column = ins_cols[j];
            for (std::uint32_t k = 0; k < num_rows; k++) { // move down the column
                if (k%2==0) {
                    block[k].push_back(ins_column[k/2]);
                    if (k != 0) {
                        if (block[k].back() != '_' && block[k-2].back() != '_' && block[k].back() != block[k-2].back()) {
                            block[k-1].push_back('!');
                        } else {
                            block[k-1].push_back(' ');
                        }
                    }        
                } 
            }
        }

           
    }

    std::uint32_t counter = 0;
    //printing
    for (auto& block: output_blocks) {
        std::cout << start_of_blocks[counter++] << std::endl;
        for (auto& row: block) {
            std::cout << row << std::endl;       
        }
        std::cout << std::endl;
    }
    
    
    
}





       
} // namespace align_reads
