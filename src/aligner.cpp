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

#define MATRIX_HEIGHT 10
#define MATRIX_WIDTH 200
#define PAD_CODE 5
#define GAP_CODE 4 
#define MIN_GOOD_OVLP_RATIO 0.95
#define EXTRACT nanosim_extract
#define OVLP_THRES 0
#define TRIALS 1000
/*
 * To Do: 
 * 1. measure coverage in each exact window, estimate coverage in each estimated window, tweak threhold/accept reject based on demands in windows.
 * 2. check how many haploytype seqs are found 
 * 3. 
 */


std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace align_reads {

// ------------------------------------

static std::uint64_t true_overlap_num = 0;
struct seq_info {
    std::uint32_t id;
    bool forward;
    std::uint32_t start;
    std::uint32_t end;
    std::string contig;
    double norm_edit_distance = -1;
    double norm_score = -1;
    std::string target_name;
    std::uint32_t dist = 0;
    std::uint32_t len = 0;
    std::uint32_t score = 0;
    std::uint32_t result_idx;
    std::uint32_t total_overhang = 0;
    double total_ratio = -1;
    std::uint32_t larger_overhang = 0;
    double larger_ratio = -1;
    std::string name;
    std::uint32_t q_left_clip;
    std::uint32_t t_left_clip;
    std::uint32_t q_right_clip;
    std::uint32_t t_right_clip;
    std::uint32_t left_extend;
    std::uint32_t right_extend;
    bool q_rev;
    std::uint32_t t_start;
    std::uint32_t q_start;
    std::uint32_t t_end;
    std::uint32_t q_end;
    std::uint32_t q_idx;
    std::uint32_t t_idx;
};	

struct comp {
    bool operator() (const seq_info& lhs, const seq_info& rhs) {
        return lhs.id < rhs.id;
    }	
};


static void print_seq_info(const seq_info& info) {
    std::cout << "id " <<  info.id << std::endl;
    std::cout << "contig " <<  info.contig << std::endl;
    std::cout << "start " <<  info.start << std::endl;
    std::cout << "end " << info.end << std::endl;
    std::cout << "forward ? " << info.forward << std::endl;
}


std::set<seq_info, comp> all_overlaps;
	
static seq_info nanosim_extract(std::string name) {
    //>NC-004354_15504699;aligned_0_R_9_1481_32
    //>NC-004354_8681245_unaligned_17149_F_0_1489_0
        
    seq_info info;
    info.name = name;
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



void Aligner::test() {
      
    std::set<align_reads::seq_info, comp> found;
    for (auto i = 0; i < TRIALS; i++) {
        auto& target = sequences[i];	
        std::uint32_t target_id = target->id;
        std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, true); 
        auto target_info = EXTRACT(target->name);       
        for (auto& o : overlaps) {
            auto& s = sequences[id_to_pos_index[o.rhs_id]];
            auto query_info = EXTRACT(s->name);
            query_info.t_idx = i;
            query_info.q_idx = id_to_pos_index[o.rhs_id];
            query_info.target_name = target_info.name;
            query_info.q_start = o.rhs_begin;
            query_info.q_end = o.rhs_end;
            query_info.t_start = o.lhs_begin;
            query_info.t_end = o.lhs_end;
            //std::cout << "F " << target_info.id << " has " << query_info.id << std::endl;
            if (overlap_len(query_info, target_info) > OVLP_THRES) found.insert(query_info);
        }       
    }
    for (auto& item:found) {
        //if (item.id != 15843) continue;
        if (all_overlaps.find(item) == all_overlaps.end()) continue;       
        all_overlaps.erase(item);    
        std::cout << "FOUND " << std::endl;
        auto target_info = EXTRACT(item.target_name);
        auto& target_s = sequences[item.t_idx];
        auto& query_s = sequences[item.q_idx];
        std::uint32_t q_start = item.start;
        std::uint32_t q_end = item.end;
        std::uint32_t t_start = target_info.start;
        std::uint32_t t_end = target_info.end;

        int q_left_clip = t_start - q_start;
        int q_right_clip = q_end - t_end;
        int t_left_clip;
        int t_right_clip;
        if (q_left_clip >= 0) {
            t_left_clip = 0;
        } else {            
            t_left_clip = - q_left_clip;
            q_left_clip = 0;
        }
        if (q_right_clip >= 0) {
            t_right_clip = 0;
        } else {
            t_right_clip = -q_right_clip;
            q_right_clip = 0;
        }
        q_left_clip += item.left_extend;
        q_right_clip += item.right_extend;
        t_left_clip += target_info.left_extend;
        t_right_clip += target_info.right_extend; 
        
        if (!item.forward) {
            query_s->ReverseAndComplement();
        }
        if (!target_info.forward) {
            target_s->ReverseAndComplement();
        }
        std::cout << "qleftc " << q_left_clip << " qrightc " << q_right_clip << std::endl;
        std::cout << "tleftc " << t_left_clip << " trightc " << t_right_clip << std::endl; 
        std::cout << "qstart " << q_start << " qend " << q_end << std::endl;
        std::cout << "tstart " << t_start << " tend " << t_end << std::endl;
        std::cout << "qostart " << item.q_start << " qoend " << item.q_end << std::endl;
       
        std::cout << "tostart " << item.t_start << " toend " << item.t_end << std::endl;
        std::cout << "qn " << item.name << std::endl;
        std::cout << "tn " << target_info.name << std::endl;

        
        std::cout << "full ql " << query_s->inflated_len << std::endl;
        std::cout << "full tl " << target_s->inflated_len << std::endl; 
        int ql = query_s->inflated_len - q_left_clip - q_right_clip;
        int tl = target_s->inflated_len - t_left_clip - t_right_clip;
        std::cout << "ql " << ql << " tl " << tl << std::endl;
        if (ql < 100 && tl < 100) {
            std::cout << "PROBLEMATIC" << std::endl;
            continue;
        } else if (ql < 100) {
            ql = tl;
        } else if (tl < 100) {
            tl = ql;
        }

        auto query = query_s->InflateData(q_left_clip, ql);
        auto target = target_s->InflateData(t_left_clip, tl);
        if (query.size() <= 100 || target.size() <= 100) {
            std::cout << "PROBLEMATIC" << std::endl;
            continue;

        }
        //std::cerr << ">" << query_s->name << std::endl << query_s->InflateData(0, -1) << std::endl;
        //std::cerr << ">" << target_s->name << std::endl << target_s->InflateData(0, -1) << std::endl;

        std::vector<std::string> overlapping_reads;
        overlapping_reads.push_back(query);
        
        EdlibAlignResult result = edlibAlign(query.c_str(), query.size(),
            	        target.c_str(), target.size(),
             	        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));        
        std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
        clips.emplace_back(0, 0);

        std::vector<EdlibAlignResult> rs = {result};
        auto alignment = pseudoMSA(overlapping_reads, target, clips, rs, false);
        alignment.print(); 
        //exit(0);
        /*
        target = query_s->InflateData(0, -1);
        query = target_s->InflateData(0, -1);
        overlapping_reads.clear();
        overlapping_reads.push_back(query);
        
        result = edlibAlign(query.c_str(), query.size(),
            	        target.c_str(), target.size(),
             	        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));        
        clips.clear();
        clips.emplace_back(0, 0);

        rs = {result};
        alignment = pseudoMSA(overlapping_reads, target, clips, rs, false);
        alignment.print(); 
        exit(0);
        
        std::cout << "----------------------prefix--------------------- " << std::endl;
        query = query_s->InflateData(0, query_s->inflated_len);
        target = target_s->InflateData(0, target_s->inflated_len);
        overlapping_reads.clear();
        overlapping_reads.push_back(query);
        
        result = edlibAlign(query.c_str(), query.size(),
            	        target.c_str(), target.size(),
             	        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));        
        clips.clear();
        clips.emplace_back(0, 0);

        rs = {result};
        alignment = pseudoMSA(overlapping_reads, target, clips, rs, false);
        alignment.print(); 

        

        if (!item.forward) query_s->ReverseAndComplement();
        if (!target_info.forward) target_s->ReverseAndComplement();

        std::cout << "----------------------reverse prefix--------------------- " << std::endl;
        target = query_s->InflateData(0, query_s->inflated_len);
        query = target_s->InflateData(0, target_s->inflated_len);
        overlapping_reads.clear();
        overlapping_reads.push_back(query);
        
        result = edlibAlign(query.c_str(), query.size(),
            	        target.c_str(), target.size(),
             	        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));        
        clips.clear();
        clips.emplace_back(0, 0);

        rs = {result};
        alignment = pseudoMSA(overlapping_reads, target, clips, rs, false);
        alignment.print(); 
        */
         

        if (!item.forward) query_s->ReverseAndComplement();
        if (!target_info.forward) target_s->ReverseAndComplement();         
    }
    for (auto& item:all_overlaps) {
        auto target_info = EXTRACT(item.target_name);
        auto& target_s = sequences[item.t_idx];
        auto& query_s = sequences[item.q_idx];
       
        if (!item.forward) query_s->ReverseAndComplement();
        if (!target_info.forward) target_s->ReverseAndComplement();
        std::uint32_t q_start = item.start;
        std::uint32_t q_end = item.end;
        std::uint32_t t_start = target_info.start;
        std::uint32_t t_end = target_info.end;

        int q_left_clip = t_start - q_start;
        int q_right_clip = q_end - t_end;
        int t_left_clip;
        int t_right_clip;
        if (q_left_clip >= 0) {
            t_left_clip = 0;
        } else {            
            t_left_clip = - q_left_clip;
            q_left_clip = 0;
        }
        if (q_right_clip >= 0) {
            t_right_clip = 0;
        } else {
            t_right_clip = -q_right_clip;
            q_right_clip = 0;
        }
        q_left_clip += item.left_extend;
        q_right_clip += item.right_extend;
        t_left_clip += target_info.left_extend;
        t_right_clip += target_info.right_extend; 

        int ql = query_s->inflated_len - q_left_clip - q_right_clip;
        int tl = target_s->inflated_len - t_left_clip - t_right_clip;
      

        std::cout << "NOT FOUND " << std::endl;
        std::cout << "qleftc " << q_left_clip << " qrightc " << q_right_clip << std::endl;
        std::cout << "tleftc " << t_left_clip << " trightc " << t_right_clip << std::endl; 
        std::cout << "qstart " << q_start << " qend " << q_end << std::endl;
        std::cout << "tstart " << t_start << " tend " << t_end << std::endl;
        std::cout << "qostart " << item.q_start << " qoend " << item.q_end << std::endl;
       
        std::cout << "tostart " << item.t_start << " toend " << item.t_end << std::endl;
        std::cout << "qn " << item.name << std::endl;
        std::cout << "tn " << target_info.name << std::endl;

        
        std::cout << "full ql " << query_s->inflated_len << std::endl;
        std::cout << "full tl " << target_s->inflated_len << std::endl; 

       
        std::cout << "ql " << ql << " tl " << tl << std::endl;
        if (ql < 100 && tl < 100) {
            std::cout << "PROBLEMATIC" << std::endl;
            continue;
        } else if (ql < 100) {
            ql = tl;
        } else if (tl < 100) {
            tl = ql;
        }
        if (ql > 800 || tl > 800) {
            std::cout << "LONG " << std::endl;
            //std::cerr << ">" << query_s->name << std::endl << query_s->InflateData(0, -1) << std::endl;
            //std::cerr << ">" << target_s->name << std::endl << target_s->InflateData(0, -1) << std::endl;
        }
        std::cout << "-----------------infix------------------- " << std::endl;
        auto query = query_s->InflateData(q_left_clip, ql);
        auto target = target_s->InflateData(t_left_clip, tl);
        if (query.size() <= 100 || target.size() <= 100) {
            std::cout << "PROBLEMATIC" << std::endl;
            continue;

        }


        std::vector<std::string> overlapping_reads;
        if (!item.forward) {
            query_s->ReverseAndComplement();

        }
        if (!target_info.forward) {

           target_s->ReverseAndComplement();
        }
        overlapping_reads.push_back(query);
        
        EdlibAlignResult result = edlibAlign(query.c_str(), query.size(),
            	        target.c_str(), target.size(),
             	        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));        
        std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
        clips.emplace_back(0, 0);

        std::vector<EdlibAlignResult> rs = {result};
        auto alignment = pseudoMSA(overlapping_reads, target, clips, rs, false);
        alignment.print(); 
        /* 

        std::cout << "----------------reverse infix---------------- " << std::endl;
        target = query_s->InflateData(0, query_s->inflated_len);
        query  = target_s->InflateData(0, target_s->inflated_len);
        overlapping_reads.clear();
        if (!item.forward) {
            query_s->ReverseAndComplement();

        }
        if (!target_info.forward) {

           target_s->ReverseAndComplement();
        }
        overlapping_reads.push_back(query);
        
        result = edlibAlign(query.c_str(), query.size(),
            	        target.c_str(), target.size(),
             	        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));        
        clips.clear();
        clips.emplace_back(0, 0);

        rs = {result};
        alignment = pseudoMSA(overlapping_reads, target, clips, rs, false);
        alignment.print(); 





        std::cout << "-----------------------prefix--------------------- " << std::endl;
        query = query_s->InflateData(0, query_s->inflated_len);
        target = target_s->InflateData(0, target_s->inflated_len);
        overlapping_reads.clear();
        if (!item.forward) {
            query_s->ReverseAndComplement();

        }
        if (!target_info.forward) {

           target_s->ReverseAndComplement();
        }
        overlapping_reads.push_back(query);
        
        result = edlibAlign(query.c_str(), query.size(),
            	        target.c_str(), target.size(),
             	        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));        
        clips.clear();
        clips.emplace_back(0, 0);

        rs = {result};
        alignment = pseudoMSA(overlapping_reads, target, clips, rs, false);
        alignment.print(); 

        std::cout << "-----------------------reverse prefix--------------------- " << std::endl;
        target = query_s->InflateData(0, query_s->inflated_len);
        query = target_s->InflateData(0, target_s->inflated_len);
        overlapping_reads.clear();
        if (!item.forward) {
            query_s->ReverseAndComplement();

        }
        if (!target_info.forward) {

           target_s->ReverseAndComplement();
        }
        overlapping_reads.push_back(query);
        
        result = edlibAlign(query.c_str(), query.size(),
            	        target.c_str(), target.size(),
             	        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));        
        clips.clear();
        clips.emplace_back(0, 0);
          
        rs = {result};
        alignment = pseudoMSA(overlapping_reads, target, clips, rs, false);
        alignment.print(); 
        */  
        if (!item.forward) {
            query_s->ReverseAndComplement();

        }
        if (!target_info.forward) {

           target_s->ReverseAndComplement();
        }


    }
   	
}

void Aligner::test2() {
    std::uint64_t tp = 0;
    std::uint64_t fp = 0;
    std::uint64_t tn = 0;
    std::uint64_t fn = 0;
    std::uint64_t total = 0;
    std::uint64_t total_overlaps = 0;
    std::uint32_t passed = 0;
    //double total_dist_tp = 0;
    //double total_dist_fn = 0;
    //double total_dist_fp = 0;
    //double total_dist_
    //double total_score_tp = 0;
    //double total_score_fp = 0;
    std::vector<align_result> results;

    for (auto i = 0; i < TRIALS; i++) {
        std::set<align_reads::seq_info, comp> all;
        std::set<align_reads::seq_info, comp> accepts;        
        auto& target = sequences[i];	
        // get target string 
        //std::string target_string = target->InflateData();      
        std::uint32_t target_id = target->id;
        //std::cout << "**************************************************************************" << std::endl;
        //std::cout << "target name: " << target->name << std::endl;
        //std::cout << "own len " << target->inflated_len << std::endl;
        // find all overlaps with target
        std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, true); 
        //std::cout << "ov num " << overlaps.size() << std::endl;        
       
        //if (overlaps.size() > 1000) continue;
        //std::cout << "ov num " << overlaps.size() << std::endl;
        //auto sort_by_id = [] (const biosoup::Overlap& o1, const biosoup::Overlap& o2) {
         //   return o1.rhs_id < o2.rhs_id;        
        //};
        
        // sort overlaps by id  
        //std::sort(overlaps.begin(), overlaps.end(), sort_by_id);
        // Fill up infos and overlapping reads 
        
        
        double prop = 0.055;
        
        auto target_info = EXTRACT(target->name); 
        for (auto& o : overlaps) {
            if (true) {
                auto& s = sequences[id_to_pos_index[o.rhs_id]];
                auto query_info = EXTRACT(s->name);
                std::uint32_t t_begin = o.lhs_begin;
                std::uint32_t t_end = o.lhs_end;
                std::uint32_t q_begin = o.rhs_begin;
                std::uint32_t q_end = o.rhs_end;

                if (!o.strand) {
                    s->ReverseAndComplement();
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

                std::uint32_t left_overhang = q_begin - q_clip_left;
                std::uint32_t right_overhang = s->inflated_len - q_clip_right - q_end;
                double bad_left_ratio = (double) left_overhang / (s->inflated_len - q_clip_left - q_clip_right);
                double bad_right_ratio = (double) right_overhang / (s->inflated_len - q_clip_left - q_clip_right);
                
                double larger_ratio = std::max(bad_left_ratio, bad_right_ratio);
	            
                bool accept = true;
                if (larger_ratio > prop) {
                    accept = false;
                } 
                std::string target_string = target->InflateData(t_clip_left, target->inflated_len - t_clip_left - t_clip_right);
                std::string query = s->InflateData(q_clip_left, s->inflated_len - q_clip_left - q_clip_right);
                EdlibAlignResult result = edlibAlign(query.c_str(), query.size(),
            	        target_string.c_str(), target_string.size(),
             	        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                
                query_info.q_left_clip = q_clip_left;
                query_info.q_right_clip = q_clip_right;
                query_info.t_left_clip = t_clip_left;
                query_info.t_right_clip = t_clip_right;
                query_info.q_rev = !o.strand;
                query_info.norm_edit_distance = (double) result.editDistance / query.size();
                query_info.norm_score = (double) o.score / query.size();
                query_info.len = s->inflated_len;
                query_info.target_name = target->name;
                query_info.dist = result.editDistance;
                query_info.score = o.score;
                std::vector<std::string> overlapping_reads = {query};
                std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
                clips.emplace_back(0, 0);
                std::vector<EdlibAlignResult> rs = {result};
                auto alignment = pseudoMSA(overlapping_reads, target_string, clips, rs, false);
                query_info.result_idx = results.size();
                results.push_back(std::move(alignment));
                query_info.q_start = q_begin;
                query_info.t_start = t_begin;
                query_info.larger_ratio = larger_ratio;
                query_info.total_ratio = bad_left_ratio + bad_right_ratio;
                query_info.larger_overhang = std::max(left_overhang, right_overhang);
                query_info.total_overhang = left_overhang + right_overhang;  
                assert(query_info.larger_overhang < 1000000);

                all.insert(query_info);
                if (accept) {
                    accepts.insert(query_info);
                } 

                if (!o.strand) s->ReverseAndComplement();
            }
        }
        double ratio = all.size() / (double) unique_counts[i];
        //if (ratio > 1.01) continue;
        passed++;
        //std::cout << "all size " << all.size() << std::endl;

        total_overlaps += overlaps.size();
        //std::cout << "num unique reported overlapping sequences: " << all.size() << std::endl;
        //std::cout << "num unique accepted: " << accepts.size() << std::endl;

        auto print_f = [&results] (std::string t, const align_reads::seq_info& info) {
             auto t_info = EXTRACT(info.target_name);
             std::cout << "----"<< t<<  "----" << std::endl;
             std::cout << "larger ratio: " << info.larger_ratio << std::endl;
             std::cout << "q full len: " << info.len << ", q_left_clip: " << info.q_left_clip << ", q_right_clip: " << info.q_right_clip << ", t_left_clip: "<< info.t_left_clip << ", t_right_clip: " << info.t_right_clip << ", qrev?: " << info.q_rev << std::endl;
             std::cout << "normalized distance: "  <<  info.norm_edit_distance << ", normalized score: " << info.norm_score  << std::endl;
             std::cout << std::endl;
             std::cout <<  "query, contig: " << info.contig << ", start: " << info.start << ", end: " << info.end  <<std::endl;
             std::cout << "query name: " << info.name << std::endl;
             std::cout << std::endl;
             std::cout << "target, " << "contig: " << t_info.contig << ", start: " << t_info.start << ", end: " << t_info.end << std::endl;
             std::cout << "target name: " << t_info.name << std::endl; 
             std::cout << "qstart: " << info.q_start << " tstart: " << info.t_start << std::endl;
             std::cout << t<< ", " << info.norm_edit_distance << ", " << info.norm_score << std::endl;
             results[info.result_idx].print();

        };
        for (auto& info : accepts) {
             std::uint32_t ovl = overlap_len(info, target_info);	    
             std::uint16_t lv = log10(ovl);
             
             //std::cout << "overlap len: " << ovl << std::endl;
             auto t_info = EXTRACT(info.target_name);
             if (ovl > OVLP_THRES) {
                 tp++;
                 print_f("tp", info);
                 //total_dist_tp += info.norm_edit_distance;
                 //total_score_tp += info.norm_score;
                 /*std::cout << "----TP----" << std::endl;
                 std::cout << "larger ratio: " << info.larger_ratio << std::endl;
                 std::cout << "q full len: " << info.len << ", q_left_clip: " << info.q_left_clip << ", q_right_clip: " << info.q_right_clip << ", t_left_clip: "<< info.t_left_clip << ", t_right_clip: " << info.t_right_clip << ", qrev?: " << info.q_rev << std::endl;
                 std::cout << "normalized distance: "  <<  info.norm_edit_distance << ", normalized score: " << info.norm_score  << std::endl;
                 std::cout << std::endl;
                 std::cout <<  "query, contig: " << info.contig << ", start: " << info.start << ", end: " << info.end  <<std::endl;
                 std::cout << "query name: " << info.name << std::endl;
                 std::cout << std::endl;
                 std::cout << "target, " << "contig: " << t_info.contig << ", start: " << t_info.start << ", end: " << t_info.end << std::endl;
                 std::cout << "target name: " << t_info.name << std::endl; 
                 std::cout << "tp, " << info.norm_edit_distance << ", " << info.norm_score << std::endl;
                 results[info.result_idx].print();*/

             } else {

                 print_f("fp", info);
                 /*std::cout << "----FP----" << std::endl;
                 std::cout << "larger ratio: " << info.larger_ratio << std::endl;
                 std::cout << "q full len: " << info.len << ", q_left_clip: " << info.q_left_clip << ", q_right_clip: " << info.q_right_clip << ", t_left_clip: "<< info.t_left_clip << ", t_right_clip: " << info.t_right_clip << ", qrev?: " << info.q_rev << std::endl;
                 std::cout << "normalized distance: "  <<  info.norm_edit_distance << ", normalized score: " << info.norm_score  << std::endl;
                 std::cout << std::endl;
                 std::cout <<  "query, contig: " << info.contig << ", start: " << info.start << ", end: " << info.end  <<std::endl;
                 std::cout << "query name: " << info.name << std::endl;
                 std::cout << std::endl;
                 std::cout << "target, " << "contig: " << t_info.contig << ", start: " << t_info.start << ", end: " << t_info.end << std::endl;
                 std::cout << "target name: " << t_info.name << std::endl;
                 std::cout << "fp, " << info.norm_edit_distance << ", " << info.norm_score << std::endl;
                 results[info.result_idx].print();
                 //total_dist_fp += info.norm_edit_distance;                 
                 //total_score_fp += info.norm_score;
*/
                 fp++;
             }

        }
       

        total += all.size();
        for (auto& info : accepts) {
            all.erase(info);
        }
        // these are rejected
        for (auto& info : all) {
             
             auto t_info = EXTRACT(info.target_name);
             std::uint32_t ovl = overlap_len(info, target_info);	    
             //std::cout << "overlap len: " << ovl << std::endl;
             if (ovl > OVLP_THRES) {
                 fn++;     

                 print_f("fn", info);
                 /*std::cout << "----FN----" << std::endl;
                 std::cout << "larger ratio: " << info.larger_ratio << std::endl;
                 std::cout << "q full len: " << info.len << ", q_left_clip: " << info.q_left_clip << ", q_right_clip: " << info.q_right_clip << ", t_left_clip: "<< info.t_left_clip << ", t_right_clip: " << info.t_right_clip << ", qrev?: " << info.q_rev << std::endl;
                 std::cout << "normalized distance: "  <<  info.norm_edit_distance << ", normalized score: " << info.norm_score  << std::endl;
                 std::cout << std::endl;
                 std::cout <<  "query, contig: " << info.contig << ", start: " << info.start << ", end: " << info.end  <<std::endl;
                 std::cout << "query name: " << info.name << std::endl;
                 std::cout << std::endl;
                 std::cout << "target, " << "contig: " << t_info.contig << ", start: " << t_info.start << ", end: " << t_info.end << std::endl;
                 std::cout << "target name: " << t_info.name << std::endl;
                
                 std::cout << "fn, " << info.norm_edit_distance << ", " << info.norm_score << std::endl;
                 results[info.result_idx].print();*/
             } else {
                 tn++;

                 print_f("tn", info);
                 /*std::cout << "----TN----" << std::endl;
                 std::cout << "larger ratio: " << info.larger_ratio << std::endl;
                 std::cout << "q full len: " << info.len << ", q_left_clip: " << info.q_left_clip << ", q_right_clip: " << info.q_right_clip << ", t_left_clip: "<< info.t_left_clip << ", t_right_clip: " << info.t_right_clip << ", qrev?: " << info.q_rev << std::endl;
                 std::cout << "normalized distance: "  <<  info.norm_edit_distance << ", normalized score: " << info.norm_score  << std::endl;
                 std::cout << std::endl;
                 std::cout <<  "query, contig: " << info.contig << ", start: " << info.start << ", end: " << info.end  <<std::endl;
                 std::cout << "query name: " << info.name << std::endl;
                 std::cout << std::endl;
                 std::cout << "target, " << "contig: " << t_info.contig << ", start: " << t_info.start << ", end: " << t_info.end << std::endl;
                 std::cout << "target name: " << t_info.name << std::endl;
                 
                 std::cout << "tn, " << info.norm_edit_distance << ", " << info.norm_score << std::endl;
                 results[info.result_idx].print();

                // exit(0);*/
             }
        }
    }
    //for (auto i = 0; i < 6; i++) {
    //	std::cout << "cat " << i << " has " << counts[i] << std::endl;
    //}
    std::cout << "passed: " << passed << std::endl;
    std::cout << "total overlaps: " << total_overlaps << std::endl;
    std::cout << "total unique reported overlapping sequences: " << total << std::endl;
    std::cout << "recall by RAM: " << (tp + fn) / (double) true_overlap_num << std::endl;
    std::cout << "tp: " << tp << std::endl;
    std::cout << "tn: " << tn << std::endl;
    std::cout << "fp: " << fp << std::endl;
    std::cout << "fn: " << fn << std::endl;
    std::cout << "accuracy: " << (tp + tn) / (double) total << std::endl;
    std::cout << "recall: " << (double) tp/ (tp + fn) << std::endl;
    std::cout << "precision: " << (double) tp / (tp + fp) << std::endl;
    //std::cout << "avg dist tp: " << total_dist_tp/tp << std::endl;
    //std::cout << "avg dist fp: " << total_dist_fp/fp << std::endl;
    //std::cout << "avg score tp: " << total_score_tp / tp << std::endl;
    //std::cout << "avg score fp: " << total_score_fp/fp << std::endl;
    std::cout << "---------------" << std::endl;	
}



void Aligner::test3() {

}



//------------------------------------------------

Data Aligner::next() {
    /*bool has_hap = !haplotypes_sequences.empty();
    
    if (has_hap) {
        
        auto result = align_overlapping_plus_haplotypes(sequences[num_processed++]);
        if (!result.valid) return Data();
        //result.alignment.print();
        return result.produce_data(has_hap, start_of_other_phase);
    } else {
        auto result = align_overlapping(sequences[num_processed++]);
        if (!result.valid) return Data();
        //result.alignment.print();
        return result.produce_data(has_hap, start_of_other_phase);
    }*/
    test3();
    return Data(); 
}

constexpr static std::uint8_t ENCODER[] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255,    
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 0, 255, 1, 255, 255,
    255, 2, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 3, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 4, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255
};
constexpr static char DECODER[] = {
    'A', 'C', 'G', 'T', '_'
};




Aligner::Aligner(const char** sequences_paths, std::shared_ptr<thread_pool::ThreadPool>& pool, std::uint8_t kmer_len, 
    std::uint8_t window_len, double freq, const char** haplotypes_paths) 
    : pool(pool), minimizer_engine(pool, kmer_len, window_len, 500, 4, 100, 500) {
    
	
    //srand (time(NULL));
    srand(422);
    biosoup::NucleicAcid::num_objects = 0;
    auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
    };    
    for (std::uint8_t i = 0; i < 2; i++) {
        const char* path = sequences_paths[i];
        if (path == nullptr) break;
        start_of_other_phase = sequences.size();
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
    
    
    if (haplotypes_paths != nullptr) {
        haplotypes_sequences.resize(2);
        haplotypes_minimizer_engines.resize(2);
    
        for (uint8_t i = 0; i < 2; i++) {
            const char* path = haplotypes_paths[i];
            std::string seq_path {path};
            if (is_suffix(seq_path, ".fasta")) {
                try { 
                    auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(seq_path);
                    auto s = p->Parse(-1);
                    auto& the_sequences = haplotypes_sequences[i];
                    the_sequences.insert(the_sequences.end(), std::make_move_iterator(s.begin()), std::make_move_iterator(s.end()));
                } catch (const std::invalid_argument& exception) {                   
                    std::cerr << exception.what() << std::endl;
                }

            } else {
                throw std::invalid_argument("[align_reads::Aligner::Aligner] Error: Invalid sequences file format.");

            }                           
        }    
    }
	
    std::uint32_t total_seq_num = 0;
    total_seq_num += sequences.size();
    
    

    for (auto& s: haplotypes_sequences) {
        total_seq_num += s.size();
    }

    
    id_to_pos_index.resize(total_seq_num);
    //std::sort(sequences.begin(), sequences.end(), long_seq_first);
    
    // ----------------------------------
    std::random_shuffle(sequences.begin(), sequences.end());
    
    //std::cout << "num seq " << sequences.size() << std::endl;    
    
    /*std::uint32_t max_len = 0;
    std::uint32_t min_len = -1;
    std::uint64_t total_len = 0;
    for (auto& s: sequences) {
        std::uint32_t len = s->inflated_len;
        if (len > max_len) max_len = len;
	if (len < min_len) min_len = len;
	total_len += len;	

    }

    std::cout << "min " << min_len << " max " << max_len << " avg " << total_len / (double) sequences.size() << std::endl;
    */

    /*
    std::cout << "full string " << sequences[0]->name << std::endl;
    auto info = EXTRACT(sequences[0]->name);
    std::cout << "id " <<  info.id << std::endl;
    std::cout << "contig " <<  info.contig << std::endl;
    std::cout << "start " <<  info.start << std::endl;
    std::cout << "end " << info.end << std::endl;
    std::cout << "forward ? " << info.forward << std::endl;
    */
    
    /*std::uint32_t ok = 0;
    for (auto& s: sequences) {
        if (ok == 1000) exit(0);
        std::cout << ">" << s->name << std::endl;
        std::cout << s->InflateData(0, -1) << std::endl;
        ok++;
    }
    exit(0);*/
    

    for (auto i = 0; i < TRIALS; i++) {        

        auto target_info = EXTRACT(sequences[i]->name);
        std::uint32_t c = 0;
        for (auto& s: sequences) {
            auto query_info = EXTRACT(s->name);
            query_info.t_idx = i;
            query_info.q_idx = c++;
            query_info.name = s->name;
            query_info.target_name = sequences[i]->name;
            if (query_info.id == target_info.id) continue;
            std::uint32_t len = overlap_len(query_info, target_info);

            if (len > OVLP_THRES) {
                true_overlap_num++;
                all_overlaps.insert(query_info);
                //std::cout << "T " << target_info.id << " has " << query_info.id << std::endl;
            }    
        } 
    }
    std::cout << "overlap num " << num_overlap << std::endl;
    


    //for (auto i = 0; i < 6; i++) {
    //	std::cout << "cat " << i << " " << counts[i] << std::endl;
    //}	    
    // -------------------------------------------

    std::uint32_t pos_index = 0;

    for (auto& s: sequences) {
        id_to_pos_index[s->id] = pos_index++; 
    }
	
    for (auto& seqs : haplotypes_sequences) {
        pos_index = 0;
        for (auto& s: seqs) {
            id_to_pos_index[s->id] = pos_index++;    
        }
    }
    
    minimizer_engine.Minimize(sequences.begin(), sequences.end(), true);
    minimizer_engine.Filter(freq);
    for (std::uint8_t i = 0; i < haplotypes_sequences.size(); i++) {
        auto& seqs = haplotypes_sequences[i];
        auto& m = haplotypes_minimizer_engines[i];
        m.Minimize(seqs.begin(), seqs.end(), true);
        m.Filter(freq);       
    }    
	
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
                   
                    // insertions in the query to the sides of target is ignored
                    if (next_target_index == 0) {
                        contained = false;
                        next_query_index++;
                        continue;
                    } else if (next_target_index == target.size()) {
                        contained = false;
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
        result.align_boundaries.emplace_back(align_start, align_end, contained);
        
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
        
        result.align_boundaries.emplace_back(align_start, align_end, contained);
        
        edlibFreeAlignResult(align); 
    }
    
    return result;
}




Aligner::align_result Aligner::pseudoMSA(std::vector<std::string>& queries, std::string& target, 
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

Aligner::align_overlapping_result Aligner::align_overlapping(std::unique_ptr<biosoup::NucleicAcid>& target) {

    // get target string 
    std::string target_string = target->InflateData();      
    std::uint32_t target_id = target->id;
    // find all overlaps with target
    std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, true); 
    auto sort_by_id = [] (const biosoup::Overlap& o1, const biosoup::Overlap& o2) {
        return o1.rhs_id < o2.rhs_id;        
    };
    
    // sort overlaps by id  
    std::sort(overlaps.begin(), overlaps.end(), sort_by_id);
    // Fill up infos and overlapping reads 
    
    
    //double prop = pa[kk];
    double prop = 0.055;
    

    std::vector<seq_info> infos;
    infos.reserve(overlaps.size());

    std::vector<std::string> overlapping_reads;
    overlapping_reads.reserve(overlaps.size());

    std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
    clips.reserve(overlaps.size());

    std::vector<EdlibAlignResult> edlib_results;
    edlib_results.reserve(overlaps.size());
    std::uint64_t last_id = -1;


    auto thread_task = [] (std::string& query, std::string& target, std::uint32_t t_clip_left, std::uint32_t t_clip_right) {
	EdlibAlignResult result = edlibAlign(query.c_str(), query.size(),
            	target.c_str() + t_clip_left, target.size() - t_clip_left - t_clip_right,
            	edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	return result;
    };

    std::vector<std::future<EdlibAlignResult>> futures;

   
    for (auto& o : overlaps) {
        if (o.rhs_id != last_id) {
            auto& s = sequences[id_to_pos_index[o.rhs_id]];
            std::uint32_t t_begin = o.lhs_begin;
            std::uint32_t t_end = o.lhs_end;
            std::uint32_t q_begin = o.rhs_begin;
            std::uint32_t q_end = o.rhs_end;
            if (!o.strand) {
                s->ReverseAndComplement();
                q_end = s->inflated_len - o.rhs_begin;
                q_begin = s->inflated_len - o.rhs_end;
            }
            int protrude_left = q_begin - t_begin + 10;
            int protrude_right = (s->inflated_len - q_end) - (target_string.size() - t_end) + 10; 

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

	    std::uint32_t left_overhang = q_begin - q_clip_left;
	    std::uint32_t right_overhang = s->inflated_len - q_clip_right - q_end;
	    double bad_left_ratio = (double) left_overhang / (s->inflated_len - q_clip_left - q_clip_right);
	    double bad_right_ratio = (double) right_overhang / (s->inflated_len - q_clip_left - q_clip_right);
	    
	    std::uint32_t larger = std::max(left_overhang, right_overhang);
	    double larger_ratio = std::max(bad_left_ratio, bad_right_ratio);
	    

	    if (larger_ratio > prop) {
	        if (!o.strand) s->ReverseAndComplement();
	        continue;
		
	    } 

            last_id = o.rhs_id;
            infos.emplace_back(o.rhs_id, !o.strand);                
            clips.emplace_back(t_clip_left, t_clip_right);
	    std::string query = s->InflateData(q_clip_left, s->inflated_len - q_clip_left - q_clip_right);
            overlapping_reads.push_back(query);
            futures.emplace_back(pool->Submit(thread_task, query, target_string, t_clip_left, t_clip_right)); 
	    if (!o.strand) s->ReverseAndComplement();
        }
    }
    std::uint32_t i = 0;
    for (auto& f: futures) {
        auto result = f.get();
	double score = result.editDistance/ (double) overlapping_reads[i++].size();		
    }

    //move behind
    for (auto r : edlib_results) {
        edlibFreeAlignResult(r);
    } 
       	
    return align_overlapping_result(); 
    //auto alignment = pseudoMSA(overlapping_reads, target_string, clips, edlib_results, false);   
   //return align_overlapping_result(std::move(alignment), std::move(infos), target_id);
}


Aligner::align_overlapping_result Aligner::align_overlapping_plus_haplotypes(std::unique_ptr<biosoup::NucleicAcid>& target) {
     // get target string 
    std::string target_string = target->InflateData();      
    std::uint32_t target_id = target->id;
    // find all overlaps with target
    std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, true); 
    auto sort_by_id = [] (const biosoup::Overlap& o1, const biosoup::Overlap& o2) {
        return o1.rhs_id < o2.rhs_id;        
    };
    // sort overlaps by id  
    std::sort(overlaps.begin(), overlaps.end(), sort_by_id);
    // Fill up infos and overlapping reads 
    std::vector<seq_info> infos;
    infos.reserve(overlaps.size());

    std::vector<std::string> overlapping_seqs;
    overlapping_seqs.reserve(overlaps.size() + 2);

    std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
    clips.reserve(overlaps.size() + 2);

    std::vector<EdlibAlignResult> edlib_results;
    edlib_results.reserve(overlaps.size());

    std::uint64_t last_id = -1;

    for (auto& o : overlaps) {
        if (o.rhs_id != last_id) {
            auto& s = sequences[id_to_pos_index[o.rhs_id]];
            std::uint32_t t_begin = o.lhs_begin;
            std::uint32_t t_end = o.lhs_end;
            std::uint32_t q_begin = o.rhs_begin;
            std::uint32_t q_end = o.rhs_end;
            if (!o.strand) {
                s->ReverseAndComplement();
                q_end = s->inflated_len - o.rhs_begin;
                q_begin = s->inflated_len - o.rhs_end;
            }
            int protrude_left = q_begin - t_begin;
            int protrude_right = (s->inflated_len - q_end) - (target_string.size() - t_end); 
            
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
            
          
            float good_overlap_ratio = (o.rhs_end - o.rhs_begin) / (float) (s->inflated_len - q_clip_left - q_clip_right);
	    if (good_overlap_ratio < MIN_GOOD_OVLP_RATIO) {
	        if (!o.strand) s->ReverseAndComplement();
		continue;
	    }
	    last_id = o.rhs_id;
            infos.emplace_back(o.rhs_id, !o.strand);
	    clips.emplace_back(t_clip_left, t_clip_right);
	    std::string query = s->InflateData(q_clip_left, s->inflated_len - q_clip_left - q_clip_right);
            overlapping_seqs.push_back(query);
            EdlibAlignResult align = edlibAlign(query.c_str(), query.size(),
            	target_string.c_str() + t_clip_left, target_string.size() - t_clip_left - t_clip_right,
            	edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	    edlib_results.push_back(align);
            if (!o.strand) s->ReverseAndComplement();

        }
    }
    
    for (std::uint8_t i = 0; i < 2; i++) {
        // find overlap with each haplotype
        auto& minimizers = haplotypes_minimizer_engines[i];
        std::vector<biosoup::Overlap> overlaps = minimizers.Map(target, true, false, true);
        if (overlaps.empty()) return align_overlapping_result();
        std::uint32_t best_index = 0;
        std::uint32_t best_score = 0;
        if (overlaps.size() > 1) {
            for (std::uint32_t j = 0; j < overlaps.size(); j++) {
                auto score = overlaps[j].score;
                if (score > best_score) {
                    best_score = score;
                    best_index = j;
                }
            }                
        }
        
        
        auto& best_match = overlaps[best_index];
        std::uint32_t best_match_id = best_match.rhs_id;
        bool best_match_strand = best_match.strand;
        
        
        
        std::uint32_t adjusted_begin = best_match.rhs_begin - best_match.lhs_begin;
        std::uint32_t adjusted_end = target_string.size() - best_match.lhs_end + best_match.rhs_end;
        
        // Extra allowance 
        //adjusted_begin = adjusted_begin < 50 ? 0 : adjusted_begin - 50;
        //adjusted_end += 50;
        
        auto& seqs = haplotypes_sequences[i];
        auto& seq = seqs[id_to_pos_index[best_match_id]];
        if (!best_match_strand) {
            seq->ReverseAndComplement();
        }
        std::string hap_string = seq->InflateData(adjusted_begin, adjusted_end - adjusted_begin);
        clips.emplace_back(0, 0);
        overlapping_seqs.push_back(hap_string);
        infos.emplace_back(best_match_id, !best_match_strand);
        EdlibAlignResult align = edlibAlign(hap_string.c_str(), hap_string.size(),
            target_string.c_str(), target_string.size(),
            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	edlib_results.push_back(align);


        if (!best_match_strand) seq->ReverseAndComplement();                
    }

    
    auto alignment = pseudoMSA(overlapping_seqs, target_string, clips, edlib_results, true);    
    
    return align_overlapping_result(std::move(alignment), std::move(infos), target_id);        
}






Data Aligner::align_overlapping_result::produce_data(bool produce_labels, std::uint32_t start_of_other_phase) {
    
    Data d;
    npy_intp dims[2];
    dims[0] = MATRIX_HEIGHT;
    dims[1] = MATRIX_WIDTH;
    
    auto& alignment = this->alignment;
    auto& target_columns = alignment.target_columns;
    auto& ins_columns = alignment.ins_columns;
    
    std::uint32_t alignment_width = alignment.width;
    
    std::vector<align_boundary>& boundaries = alignment.align_boundaries; 
    
    // to store window to reads info
    std::uint32_t num_windows = ceil((float) alignment_width/MATRIX_WIDTH);
    std::vector<std::vector<std::uint32_t>> window_to_reads; 
    window_to_reads.resize(num_windows);
    
    // map window to first index in it
    std::vector<std::pair<std::uint32_t, std::uint32_t>> window_to_starts;
    window_to_starts.reserve(num_windows);
    
    std::uint32_t target_index = 0;
    std::uint32_t ins_index = 0;
    std::uint32_t width_index = 0;
    std::uint32_t current_window = 0;
    
    while (width_index < alignment_width) {
    // map window to first_index 
        
        // invariant: target_index, ins_index, width_index updated at this point
        
        //calculate window 
        std::uint32_t window = width_index / MATRIX_WIDTH;               
        bool new_window_or_first_step = window > current_window || width_index == 0;
        // if 1st window or 1st step, update record first index for window
        if (new_window_or_first_step) {
            if (window > current_window) current_window++;
            window_to_starts.emplace_back(target_index, ins_index);            
        }
        
        // go through reads
        for (std::uint32_t b_i = 0; b_i < boundaries.size(); b_i++) {
            auto& b = boundaries[b_i];
            
            // check if starting point in 1st window, record read to windows and window to reads
            if (target_index == b.align_start && ins_index == 0) {
                window_to_reads[current_window].push_back(b_i);                   
            }
            
            // record for start of window after 1st(which is taken care of by above)
            if(new_window_or_first_step && target_index > b.align_start && target_index <= b.align_end) {
                window_to_reads[window].push_back(b_i);                                               
            }
           
        }
        width_index++;
        if (ins_index < ins_columns[target_index].size()) {
            ins_index++;
            
        } else {
            target_index++;
            ins_index = 0;
        }            
    }
    std::cout << "window coverage " << std::endl;
    for (auto i = 0; i < window_to_reads.size(); i++) {
	    auto& p = window_to_starts[i];
	    std::cout << " window i starts at " << p.first << ", " << p.second << std::endl;
	    std::cout << " window " << i  << " has " << window_to_reads[i].size() << " reads " << std::endl;

    }	    

    auto fill_X_from_target = [&window_to_starts, &target_columns, &ins_columns] (
        PyObject* matrix, std::uint32_t window_index, std::uint32_t row_index) {
        
        uint8_t* value_ptr;
        std::uint16_t col_index = 0;
        auto& start_index = window_to_starts[window_index];
        std::uint32_t target_index = start_index.first;
        std::uint32_t ins_index = start_index.second;
        char base;    
         
        for ( ; col_index < MATRIX_WIDTH; col_index++) {
            if (target_index < target_columns.size()) {
                if (ins_index == 0) {
                    base = target_columns[target_index][0];
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];    
                    
                } else {
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = GAP_CODE;                      
                }              
            } else {
                value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                *value_ptr = PAD_CODE;
                                
            }
                       
            if (ins_index < ins_columns[target_index].size()) {
                ins_index++;         
            } else {
                target_index++;
                ins_index = 0;
            }                            
        }   
    };
    
    auto fill_X_from_read = [&boundaries, &target_columns, &ins_columns, &window_to_starts] (PyObject* matrix, std::uint32_t window_index,
        std::uint32_t row_index, std::uint32_t read_index) {
        uint8_t* value_ptr;
        std::uint16_t col_index = 0;
        auto& start_index = window_to_starts[window_index];
        std::uint32_t target_index = start_index.first;
        std::uint32_t ins_index = start_index.second;
        char base;
              
        // fill from non target read
            auto& b = boundaries[read_index];
            for ( ; col_index < MATRIX_WIDTH; col_index++) {
                if (target_index >= b.align_start && target_index <= b.align_end) {
                    if (ins_index == 0) {
                        base = target_columns[target_index][read_index+1];
                        value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                        *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];    
                    } else {
                        base = ins_columns[target_index][ins_index - 1][read_index + 1];
                        value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                        *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];                                        
                    }                   
                                              
                } else {
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = PAD_CODE;                
                }
                if (ins_index < ins_columns[target_index].size()) {
                    ins_index++;         
                } else {
                    target_index++;
                    ins_index = 0;
                }     
            }
        
    };
    
    auto fill_Y_for_target = [&window_to_starts, &target_columns, &ins_columns, this, start_of_other_phase] (
        PyObject* matrix, std::uint32_t window_index, std::uint32_t row_index) {
        
        std::uint8_t hap = this->target_id < start_of_other_phase ? 0 : 1;    
        uint8_t* value_ptr;
        std::uint16_t col_index = 0;
        auto& start_index = window_to_starts[window_index];
        std::uint32_t target_index = start_index.first;
        std::uint32_t ins_index = start_index.second;
        char base;    
        for ( ; col_index < MATRIX_WIDTH; col_index++) {
            if (target_index < target_columns.size()) {
                if (ins_index == 0) {
                    base = target_columns[target_index][this->infos.size() - 2 + hap + 1];
                                     
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];    
                    
                } else {
                    base = ins_columns[target_index][ins_index - 1][this->infos.size() - 2 + hap + 1];
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];                      
                }              
            } else {
                value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                *value_ptr = PAD_CODE;
                                
            }
                       
            if (ins_index < ins_columns[target_index].size()) {
                ins_index++;         
            } else {
                target_index++;
                ins_index = 0;
            }                            
        }   
    };
    
    auto fill_Y_for_read = [&window_to_starts, &boundaries, &target_columns, &ins_columns, this, start_of_other_phase] (PyObject* matrix,
        std::uint32_t window_index, std::uint32_t row_index, std::uint32_t read_index) {
        uint8_t* value_ptr;
        std::uint16_t col_index = 0;
        auto& start_index = window_to_starts[window_index];
        std::uint32_t target_index = start_index.first;
        std::uint32_t ins_index = start_index.second;
        char base;
        std::uint8_t hap = this->infos[read_index].id < start_of_other_phase ? 0 : 1;       
        // fill from non target read
            auto& b = boundaries[read_index];
            for ( ; col_index < MATRIX_WIDTH; col_index++) {
                if (target_index >= b.align_start && target_index <= b.align_end) {
                    if (ins_index == 0) {
                        base = target_columns[target_index][boundaries.size() - 2 + hap +1];
                        value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                        *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];    
                    } else {
                        base = ins_columns[target_index][ins_index - 1][boundaries.size() - 2 + hap +1];
                        value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                        *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];                                        
                    }                   
                                              
                } else {
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = PAD_CODE;                
                }
                if (ins_index < ins_columns[target_index].size()) {
                    ins_index++;         
                } else {
                    target_index++;
                    ins_index = 0;
                }     
            }
        
    };

    std::uint8_t offset = produce_labels ? 2 : 0; // to exclude haplotype sequences
    
    for (std::uint32_t window_index = 0; window_index < num_windows; window_index++) {
        
        // find the reads that can be placed into this window/matrix
        std::vector<std::uint32_t>& eligible_reads = window_to_reads[window_index];
        std::sort(eligible_reads.begin(), eligible_reads.end());
        std::uint16_t num_choices = eligible_reads.size() < offset ? 0 : eligible_reads.size() - offset;
        if (num_choices == 0) continue;
	    
        // fill 1st row of matrix with target
        auto x = PyArray_SimpleNew(2, dims, NPY_UINT8);        
		fill_X_from_target(x, window_index, 0);
        
        PyObject* y = nullptr;
        if (produce_labels) {
            y = PyArray_SimpleNew(2, dims, NPY_UINT8);
            fill_Y_for_target(y, window_index, 0);
        }
        // fill the others with sampled reads
        for (std::uint32_t row_index = 1; row_index < MATRIX_HEIGHT; row_index++) {           
            auto randomn = rand() % num_choices;
            fill_X_from_read(x, window_index, row_index, eligible_reads[randomn]);
            if (produce_labels) fill_Y_for_read(y, window_index, row_index, eligible_reads[randomn]);
        }
        
        // output the x matrix    
        d.X.push_back(x);
        if (produce_labels) d.Y.push_back(y);
    }
    
    
       
    return d;
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
