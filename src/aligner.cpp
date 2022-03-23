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

#include "align_reads/aligner.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "edlib.h"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

#define MATRIX_HEIGHT 5
#define MATRIX_WIDTH 20
#define PAD_CODE 5
#define GAP_CODE 4 

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace align_reads {

Data Aligner::next() {
    auto result = align_overlapping_plus_haplotypes(sequences[num_processed++]);
    // !!! why only 20 windows? 
    // !!! keep unknown in X as unknown
    return result.produce_data(start_of_other_phase != static_cast<std::uint32_t>(-1), start_of_other_phase);
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
   
Aligner::Aligner(const char** sequences_paths, std::uint32_t num_threads, std::uint8_t kmer_len, 
    std::uint8_t window_len, double freq, const char** haplotypes_paths) 
    : minimizer_engine(std::make_shared<thread_pool::ThreadPool>(num_threads), kmer_len, window_len) {
    
	
    srand (time(NULL));
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
    
	std::sort(sequences.begin(), sequences.end(), long_seq_first);
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

Aligner::align_result Aligner::align_to_target(std::vector<std::string>& queries, std::string& target,
    bool clip_query, std::vector<std::pair<std::uint32_t, std::uint32_t>>* pads) {
    
    EdlibEqualityPair additionalEqualities[4] = {{'A', 'X'}, {'C', 'X'}, {'G', 'X'}, {'T', 'X'}};     
    auto mode = clip_query ? EDLIB_MODE_HW : EDLIB_MODE_NW;
    align_result result;
    std::uint8_t for_front_ins = clip_query ? 0 : 1;
    // Fill in the bases from the target into the first row
    result.target_columns.resize(target.size());        
    for (std::uint32_t i = 0; i < result.target_columns.size(); i++) { // for each column
        auto& column = result.target_columns[i];
        column.resize(queries.size() + 1, '_');
        column[0] = target[i]; // the first item of each row - makes up the first row
    }

    // Allowance for insertion columns after each target position
    result.ins_columns.resize(target.size() + for_front_ins);
    
    // Initialize width of the alignment to the length of the target 
    result.width = target.size(); // will be incremented along with adding insertion columns
    
    // Initialize the record of who is inserting at each position
    result.inserters.resize(target.size() + for_front_ins);
 
    
    // Aligning queries to target and filling up/adding columns
    for (std::uint32_t k = 0; k < queries.size(); k++) {      
    
        const char* padded_target_input = target.c_str();
        std::uint32_t padded_target_size = target.size();
        std::uint32_t left_pad = 0;
        std::uint32_t right_pad = 0;
        std::string padded_target;
        std::uint8_t eq_count = 0;
        if (pads != nullptr && !pads->empty()) {
            // padding
            std::pair<std::uint32_t, std::uint32_t>& p = pads->at(k); 
            left_pad = p.first;
            right_pad = p.second;
            
            
            if (left_pad + right_pad > 0) {
                padded_target.clear();
                padded_target.reserve(left_pad + target.size() + right_pad);
                padded_target.append(left_pad, 'X');
                padded_target.append(target);
                padded_target.append(right_pad, 'X');    
                padded_target_input = padded_target.c_str();
                padded_target_size = padded_target.size();
                eq_count = 4;
            } else {
                padded_target_input = target.c_str();
                padded_target_size = target.size();
                eq_count = 0;
            }                
        }
        

        auto& query = queries[k];
        EdlibAlignResult align = edlibAlign(query.c_str(), query.size(),
            padded_target_input, padded_target_size, edlibNewAlignConfig(-1, mode, EDLIB_TASK_PATH, additionalEqualities, eq_count));
        std::uint32_t t_pointer = align.startLocations[0];
        std::uint32_t next_query_index = 0; // upcoming query position
        std::uint32_t next_target_index = align.startLocations[0] <= left_pad ? 0 :align.startLocations[0] - left_pad; // upcoming target position
        std::uint32_t next_ins_index = 0; // upcoming ins index for the current segment
        std::uint32_t align_start = next_target_index;
        std::uint32_t align_end = static_cast<std::uint32_t>(align.endLocations[0] >= left_pad + target.size() ? target.size() - 1 : align.endLocations[0] - left_pad );
        bool contained = true;
        
        for (int i = 0; i < align.alignmentLength; i++) {
            
         
            switch (align.alignment[i]) {
                case 0: { // match
                     // if in padded area
                    if (t_pointer < left_pad) {
                        contained = false;
                        t_pointer++;
                        next_query_index++;
                        continue;
                    } else if (t_pointer >= left_pad + target.size()) {
                        contained = false;
                        i = align.alignmentLength;  
                        continue;   
                    }
 
                    //place matching query base from the query into column
                    
                    result.target_columns[next_target_index++][k+1] = query[next_query_index++]; 

                    next_ins_index = 0;
                    t_pointer++;
                    break;
                }
                case 1: { // ins
                   
                    std::uint32_t ins_columns_index = next_target_index - 1;
                    if (clip_query) {
                        // insertions in the query to the sides of target is ignored
                        if (next_target_index == 0) {
                            contained = false;
                            next_query_index++;
                            continue;
                        } else if (next_target_index == left_pad + target.size()) {
                            contained = false;
                            i = align.alignmentLength;  // if at ins of the query to the right of target, terminate loop early
                            continue;                       
                        }    
                    } else {
                        //Ins to left of target. 
                        if (next_target_index == 0) {
                            contained = false;
                            ins_columns_index = target.size();        
                        }
                    }
                    // Record existence of ins at certain positions
                    if (next_ins_index == 0) {
                        auto& v = result.inserters[ins_columns_index]; 
                        v.push_back(k);
                        if (v.size() == 2) {
                            result.ins_at_least2.push_back(ins_columns_index);
                        }
                    }                        
                    // check if we have enough columns ready for this target position
                    auto& ins_columns = result.ins_columns[ins_columns_index];
                    if (ins_columns.size() < next_ins_index + 1) {
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
                    t_pointer++;
                    next_ins_index = 0;
                    break;
                }
                case 3: { // mismatch
                   
                    // if in padded area
                    if (t_pointer < left_pad) {
                        contained = false;
                        t_pointer++;
                        next_query_index++;
                        continue;
                    } else if (t_pointer >= left_pad + target.size()) {
                        contained = false;
                        i = align.alignmentLength;  
                        continue;   
                    }
                    // place mismatching query base into column
                    result.target_columns[next_target_index++][k+1] = query[next_query_index++]; 
                    next_ins_index = 0;
                    t_pointer++;
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
    std::vector<std::pair<std::uint32_t, std::uint32_t>>& pads) {
    align_result result = align_to_target(queries, target, true, &pads);
    for (auto& ins_pos : result.ins_at_least2) { // for all positions that have at least two queries with insertions, 
        // extract the ins segments there
        std::uint32_t position_index_longest = 0; // index in ins_segments
        std::uint32_t longest_ins_len = 0;
        
        std::vector<std::string> ins_segments;
        auto inserters_at_pos = result.inserters[ins_pos];
        ins_segments.resize(inserters_at_pos.size());
        auto& ins_columns = result.ins_columns[ins_pos];
        for (auto& column : ins_columns) { // for each column, extract the non gap bases for the queries with ins there
            for (std::uint32_t i = 0 ; i < inserters_at_pos.size(); i++) {
                char& base = column[inserters_at_pos[i]+1];
                auto& segment = ins_segments[i];
                if (base != '_') {
                    segment.push_back(base);
                    if (segment.size() > longest_ins_len) {
                        position_index_longest = i;
                        longest_ins_len = segment.size();
                    }
                }
            }            
        }
        std::uint32_t query_id_longest_ins = inserters_at_pos[position_index_longest];
        inserters_at_pos.erase(inserters_at_pos.begin() + position_index_longest);
        std::string longest_ins_segment = std::move(ins_segments[position_index_longest]);
        ins_segments.erase(ins_segments.begin() + position_index_longest);
        auto sub_result = align_to_target(ins_segments, longest_ins_segment, false);
        
        ins_columns.resize(sub_result.width);
        std::uint32_t main_column_id = 0;
        for(std::uint32_t i = 0; i < sub_result.target_columns.size(); i++) {
            auto& ins_column = ins_columns[main_column_id++];
            auto& sub_column = sub_result.target_columns[i];
            ins_column[query_id_longest_ins + 1] = sub_column[0];
            for (std::uint32_t j = 0; j < inserters_at_pos.size(); j++) {
                ins_column[inserters_at_pos[j] + 1] = sub_column[j + 1];               
            }
            auto& sub_ins_columns = sub_result.ins_columns[i];
            for (std::uint32_t j = 0; j < sub_ins_columns.size(); j++) {
                auto& ins_column = ins_columns[main_column_id++];
                auto& sub_column = sub_ins_columns[j];
                ins_column[query_id_longest_ins + 1] = sub_column[0];    
                for (std::uint32_t j = 0; j < inserters_at_pos.size(); j++) {
                ins_column[inserters_at_pos[j] + 1] = sub_column[j + 1];               
                }    
            }            
        }           
    }

    
    return result;
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
    std::vector<std::pair<std::uint32_t, std::uint32_t>> pads;
    pads.reserve(overlaps.size() + 2);
    std::uint64_t last_id = -1;

    for (auto& o : overlaps) {
        if (o.rhs_id != last_id) {
            last_id = o.rhs_id;
            infos.emplace_back(o.rhs_id, !o.strand);
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
            std::uint32_t left_pad = q_begin > t_begin ? q_begin - t_begin : 0;
            std::uint32_t right_pad = (s->inflated_len - q_end) > (target_string.size() - t_end) ? (s->inflated_len - q_end) - (target_string.size() - t_end) : 0;  
            pads.emplace_back(left_pad, right_pad);
            overlapping_seqs.push_back(s->InflateData());
            s->ReverseAndComplement();

        }
    }
    
    for (std::uint8_t i = 0; i < 2; i++) {
        // find overlap with each haplotype
        auto& minimizers = haplotypes_minimizer_engines[i];
        std::vector<biosoup::Overlap> overlaps = minimizers.Map(target, true, false, true);
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
        adjusted_begin = adjusted_begin < 50 ? 0 : adjusted_begin - 50;
        adjusted_end += 50;
        
        auto& seqs = haplotypes_sequences[i];
        auto& seq = seqs[id_to_pos_index[best_match_id]];
        //info.emplace_back(best_match_id, !best_match_strand);
        if (!best_match_strand) {
            seq->ReverseAndComplement();
        }
        std::string hap_string = seq->InflateData(adjusted_begin, adjusted_end - adjusted_begin);
        pads.emplace_back(50, 50);
        overlapping_seqs.push_back(std::move(hap_string));
        seq->ReverseAndComplement();                
    }
    
    auto alignment = pseudoMSA(overlapping_seqs, target_string, pads);    
    return align_overlapping_result(std::move(alignment), std::move(infos), target_id);        
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
    std::vector<seq_info> infos;
    infos.reserve(overlaps.size());
    std::vector<std::string> overlapping_reads;
    overlapping_reads.reserve(overlaps.size());
    std::vector<std::pair<std::uint32_t, std::uint32_t>> pads;
    pads.reserve(overlaps.size());
    std::uint64_t last_id = -1;

    for (auto& o : overlaps) {
        if (o.rhs_id != last_id) {
            last_id = o.rhs_id;
            infos.emplace_back(o.rhs_id, !o.strand);
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
            std::uint32_t left_pad = q_begin > t_begin ? q_begin - t_begin : 0;
            std::uint32_t right_pad = (s->inflated_len - q_end) > (target_string.size() - t_end) ? (s->inflated_len - q_end) - (target_string.size() - t_end) : 0;  
            pads.emplace_back(left_pad, right_pad);
            overlapping_reads.push_back(s->InflateData());
            s->ReverseAndComplement();

        }
    }
    
    auto alignment = pseudoMSA(overlapping_reads, target_string, pads);    
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
    
    // to store read to windows info
    
    /*struct windows_info {
        std::uint32_t first_window;
        std::uint32_t last_window;
    };
        
    std::vector<windows_info> read_to_windows;
    read_to_windows.resize(boundaries.size());
    */
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
                //auto& info = read_to_windows[b_i];
                //info.first_window = current_window;
                window_to_reads[current_window].push_back(b_i);                   
            }
            
            // record for start of window after 1st(which is taken care of by above)
            if(new_window_or_first_step && target_index > b.align_start && target_index <= b.align_end) {
                window_to_reads[window].push_back(b_i);                                               
            }
            
            // check if end point in last window, record only for reads to window, window to reads taken care of by above
            //if (target_index == b.align_end && ins_index == 0) {
            //    auto& info = read_to_windows[b_i];
            //    info.last_window = current_window;    
            //}            
        }
        width_index++;
        if (ins_index < ins_columns[target_index].size()) {
            ins_index++;
            
        } else {
            target_index++;
            ins_index = 0;
        }            
    }
    
    
    /*
    // Switch to iterating through the target positions! replace below
    
    for (std::uint32_t k = 0; k < boundaries.size(); k++) {
        auto& b = boundaries[k];
        //realigning boundaries to width index
        std::uint32_t target_align_start = b.align_start; // inclusive
        std::uint32_t target_align_end = b.align_end;     // inclusive
        std::uint32_t total_ins = 0;
        for (std::uint32_t i = 0; i < target_align_start; i++) {
            auto to_add = this->alignment.ins_columns[i].size();
            // realign boundaries based on insertions
            if (to_add > 0) {
                b.align_start += to_add;
                b.align_end += to_add;
                total_ins += to_add;
            }
        }

        // where this read start in each window
        auto& starting_indices = windows_covered[k].starting_indices;
        starting_indices.reserve((b.align_end - b.align_start + 1) / MATRIX_WIDTH + 1); // estimated num of windows covered

        std::uint32_t last_window = b.align_start / MATRIX_WIDTH;
        for (std::uint32_t i = target_align_start; i <= target_align_end; i++) {
            // realign boundaries based on insertions
            b.align_end += this->alignment.ins_columns[i].size();           
            
            // location info
            std::uint32_t width_idx = i + total_ins; // where am I at on the alignment
            std::uint32_t current_window =  width_idx / MATRIX_WIDTH;
            //std::cout << "curr " << current_window << std::endl;           
            
            // if first or new window at a target position
            if (current_window > last_window || i == target_align_start) {
                // update of windows_info
                if (i == target_align_start) windows_covered[k].first_window = current_window;  // if first window              
                starting_indices.emplace_back(i, 0);
                
                // update of covering_reads
                covering_reads[current_window].push_back(k);
                
                if (i != target_align_start) last_window++;    
            }
                                                                                      
            for (std::uint32_t j = 1; j <= this->alignment.ins_columns[i].size(); j++) {
                total_ins++;
                
                // location info
                std::uint32_t width_idx = i + total_ins; // where am I at on the alignment
                std::uint32_t current_window =  width_idx / MATRIX_WIDTH;
                //std::cout << "curr " << current_window << std::endl;    
                
                // if new window at an ins position
                if (current_window > last_window) {
                    // update of windows_info                                
                    starting_indices.emplace_back(i, j);
                    
                    // update of covering_reads
                    covering_reads[current_window].push_back(k);
                    
                    last_window++;    
                }
                
            }                 
        }
        
        windows_covered[k].last_window = last_window;

    }
    */
    

    
    auto fill_row = [&window_to_reads, &window_to_starts, &boundaries, &target_columns, &ins_columns] (PyObject* matrix, std::uint32_t window_index, std::uint32_t row_index, 
        std::uint32_t read_index, bool pad, bool target) {
        uint8_t* value_ptr;
        std::uint16_t col_index = 0;
        auto& start_index = window_to_starts[window_index];
        std::uint32_t target_index = start_index.first;
        std::uint32_t ins_index = start_index.second;
        char base;
        // fill from target seq
        if (target) {
            for ( ; col_index < MATRIX_WIDTH; col_index++) {
                if (target_index >= target_columns.size()) {
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = PAD_CODE;   
                } else if (ins_index == 0) {
                    base = target_columns[target_index][0];
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];    
                } else {
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = GAP_CODE;                                        
                }
                
                if (ins_index < ins_columns[target_index].size()) {
                    ins_index++;         
                } else {
                    target_index++;
                    ins_index = 0;
                }                            
            }
            return;   
        }
        // fill with pad  
        if (pad) {            
            for (std::uint32_t col_index =0; col_index < MATRIX_WIDTH; col_index++) {
                value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                *value_ptr = PAD_CODE;
            }                    
        } else {
            // fill from non target read
            auto& b = boundaries[read_index];
            for ( ; col_index < MATRIX_WIDTH; col_index++) {
                if (target_index < b.align_start || target_index > b.align_end) {
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                    *value_ptr = PAD_CODE;                           
                } else {
                    if (ins_index == 0) {
                        base = target_columns[target_index][read_index+1];
                        value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                        *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];    
                    } else {
                        base = ins_columns[target_index][ins_index - 1][read_index + 1];
                        value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index);
                        *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];                                        
                    }               
                }
                if (ins_index < ins_columns[target_index].size()) {
                    ins_index++;         
                } else {
                    target_index++;
                    ins_index = 0;
                }     
            }
            
            /*
            auto& starting_index = info.starting_indices[window_index - info.first_window];
            std::uint16_t start_in_window = 0;
			std::uint16_t end_in_window = MATRIX_WIDTH; // exclusive
			
            if (window_index == info.first_window) {
				start_in_window = boundaries[read_index].align_start % MATRIX_WIDTH;
			} 
            
			if (window_index == info.last_window) {
				end_in_window = boundaries[read_index].align_end % MATRIX_WIDTH;
			}
            std::uint16_t col_index = start_in_window;
            
            
            char base;
            //pad the sides if needed
            
            for (std::uint16_t i = 0; i < start_in_window; i++) {
                value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, i);
                *value_ptr = PAD;                
            }
            for (std::uint16_t i = end_in_window; i < MATRIX_WIDTH; i++) {
                value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, i);
                *value_ptr = PAD;                
            }
            
            
            // if start in ins segment, finish the ins and start from next main index            
            if (starting_index.second > 0) {
                for (std::uint32_t i = starting_index.second; i <= this->alignment.ins_columns[starting_index.first].size(); i++) {
                    if (col_index > end_in_window) return;
                    base = this->alignment.ins_columns[starting_index.first][i-1][read_index + 1];
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index++);
                    *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];

                }                    
                
            }
          
            for (std::uint32_t target_index = starting_index.second > 0 ? starting_index.first + 1: starting_index.first; ; target_index++) {
                if (col_index > end_in_window) return;
                base = this->alignment.target_columns[target_index][read_index + 1];
                value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index++);
                *value_ptr = ENCODER[static_cast<std::uint8_t>(base)];
                for (std::uint32_t i = 1; i <= this->alignment.ins_columns[target_index].size(); i++) {
                    if (col_index > end_in_window) return;
                    base = this->alignment.ins_columns[target_index][i-1][read_index + 1];
                    value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row_index, col_index++);
                    *value_ptr = ENCODER[static_cast<std::uint8_t>(base)]; 
                }               
            }
            */
        }       
    };
    
    std::uint8_t offset = produce_labels ? 2 : 0; // to exclude haplotype sequences
    std::uint8_t target_hap = this->target_id < start_of_other_phase ? 0 : 1;
    
    // put read segments into windows, padd unfilled for now
    for (std::uint32_t window_index = 0; window_index < num_windows; window_index++) {
        
        // find the reads that can be placed into this window/matrix
        std::vector<std::uint32_t>& eligible_reads = window_to_reads[window_index];
        std::uint16_t num_choices = eligible_reads.size() < offset ? 0 : eligible_reads.size() - offset;
        if (num_choices == 0) continue;
        
        // to record the chosen reads' ids
        std::vector<std::uint32_t> chosen_ids;
        if (produce_labels) chosen_ids.reserve(MATRIX_HEIGHT - 1);
	    
        // fill 1st row of matrix with target
        auto x = PyArray_SimpleNew(2, dims, NPY_UINT8);        
		fill_row(x, window_index, 0, 0, false, true);
        
        // fill the others with sampled reads
        for (std::uint32_t row_index = 1; row_index < MATRIX_HEIGHT; row_index++) {           
            auto randomn = rand() % num_choices;
            fill_row(x, window_index, row_index, eligible_reads[randomn], false, false);
            if (produce_labels) chosen_ids.push_back(this->infos[randomn].id);
        }
        
        // output the x matrix    
        d.X.push_back(x);
        
        
        if (produce_labels) {
            // fill with hap sequences
            auto y = PyArray_SimpleNew(2, dims, NPY_UINT8);
            fill_row(y, window_index, 0, boundaries.size() - 2  + target_hap, false, false);    
            
            for (std::uint32_t row_index = 1; row_index < MATRIX_HEIGHT; row_index++) {           
                std::uint8_t hap = chosen_ids[row_index - 1] < start_of_other_phase ? 0 : 1;
                fill_row(y, window_index, row_index, boundaries.size() - 2 + hap, false, false);
                
            }    
            d.Y.push_back(y);   
        }
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
    std::uint32_t col_index = 0;
    if (this->ins_columns.size() > this->target_columns.size()) {
        
        auto& ins_cols = this->ins_columns[this->target_columns.size()];
        for (std::uint32_t j = 0; j < ins_cols.size() ; j++) { // for each ins column here
            auto& ins_column = ins_cols[j];
            auto& block = output_blocks[col_index++/block_width];
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
        auto& block = output_blocks[col_index++/block_width];
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
            auto& block = output_blocks[col_index++/block_width];
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
        std::cout << counter << std::endl;
        counter += block_width;
        for (auto& row: block) {
            std::cout << row << std::endl;       
        }
        std::cout << std::endl;
    }
    
    
    
}





       
} // namespace align_reads