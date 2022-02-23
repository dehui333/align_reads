#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

#include "align_reads/feature_generator.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "edlib.h"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};


namespace align_reads {

constexpr static std::uint16_t MATRIX_COL_NUM = 256;
constexpr static std::uint16_t MATRIX_ROW_NUM = 256;

    
FeatureGenerator::FeatureGenerator(const char* sequences_path, std::uint32_t num_threads, std::uint8_t kmer_len, 
    std::uint8_t window_len, double freq) 
    : minimizer_engine(std::make_shared<thread_pool::ThreadPool>(num_threads), kmer_len, window_len) {
    std::string seq_path {sequences_path};
    
    auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
    };
    
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
        throw std::invalid_argument("[align_reads::FeatureGenerator::FeatureGenerator] Error: Invalid sequences file format.");
    }
    auto long_seq_first = [] (const std::unique_ptr<biosoup::NucleicAcid>& s1, const std::unique_ptr<biosoup::NucleicAcid>& s2) {
        return s1->inflated_len > s2->inflated_len;        
    };
    
    id_to_pos_index.resize(sequences.size());
    std::sort(sequences.begin(), sequences.end(), long_seq_first);
    std::uint32_t pos_index = 0;
    for (auto& s: sequences) {
        id_to_pos_index[s->id] = pos_index++; 
    }
    minimizer_engine.Minimize(sequences.begin(), sequences.end(), true);
    minimizer_engine.Filter(freq);
}

FeatureGenerator::align_result FeatureGenerator::align_to_target(std::vector<std::string> queries, std::string target) {
    align_result r;
    return r;
}

FeatureGenerator::align_result FeatureGenerator::pseudoMSA(std::vector<std::string> queries, std::string target) {
    align_result r;
    return r;
}

FeatureGenerator::align_overlapping_result FeatureGenerator::align_overlapping(std::unique_ptr<biosoup::NucleicAcid>& target) {
    /*
    // find all overlaps with target
    std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, true);  
    auto sort_by_id = [] (const biosoup::Overlap& o1, const biosoup::Overlap& o2) {
        return o1.rhs_id < o2.rhs_id;        
    };
    
    // sort overlaps by id
    std::vector<biosoup::Overlap*> unique_overlaps;  
    std::sort(overlaps.begin(), overlaps.end(), sort_by_id);
    
    // remove duplicates
    std::uint64_t last_id = -1;
    for (auto& o : overlaps) {
        if (o.rhs_id != last_id) {
            unique_overlaps.push_back(&o);
            last_id = o.rhs_id;
        }
    }
    // get target string 
    std::string target_string = target->InflateData();          

    // Initialize results
    align_result result;
    std::vector<read_info> infos;
    // infos
    infos.resize(unique_overlaps.size() + 1);
    infos[0].id = target->id;
    infos[0].same_strand = true;
    infos[0].align_start = 0;
    infos[0].align_end = target->inflated_len;
    
    // result.target_positions_pileup
    result.target_positions_pileup.resize(target_string.size());        
    for (std::uint32_t i = 0; i < result.target_positions_pileup.size(); i++) {
        auto& column = result.target_positions_pileup[i];
        column.resize(unique_overlaps.size() + 1, ENCODER['_']);
        column[0] = ENCODER[static_cast<std::uint8_t>(target_string[i])];
    }
    
    // result.ins_positions_pileup
    result.ins_positions_pileup.resize(target_string.size());
    
    // result.width
    result.width = target_string.size(); // will be incremented
    
    // Aligning with queries and filling up results
    for (std::uint32_t k = 0; k < unique_overlaps.size(); k++) {
        auto& o = unique_overlaps[k];
        std::uint32_t pos_index = id_to_pos_index[o->rhs_id];
        auto& query_seq = sequences[pos_index];
        if (!o->strand) {
            query_seq->ReverseAndComplement();
        }
        std::string query_string = query_seq->InflateData();
        EdlibAlignResult align = edlibAlign(query_string.c_str(), query_string.size(),
            target_string.c_str(), target_string.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        std::uint32_t query_index = 0; // next query position
        std::uint32_t target_index = 0; // next target position
        std::uint32_t ins_index = 0; // next ins index for the current segment
        
        for (int i = 0; i < align.alignmentLength; i++) {
            switch (align.alignment[i]) {
                case 0:  // match
                    result.target_positions_pileup[target_index++][k+1] = ENCODER[static_cast<std::uint8_t>(query_string[query_index++])];
                    break;
                case 1:  // ins
                    if (target_index == 0) {
                        continue;
                    }
                    break;
                case 2:  // del
                    break;
                case 3:  // mismatch
                    break;
                default:
                    std::cerr << "Unknown alignment result by edlib!" << std::endl;               
            }            
        }
        
        edlibFreeAlignResult(align); 
    }
    
    
    
    
    return result;*/
    align_overlapping_result r;
    return r;    
}

void FeatureGenerator::print_align(FeatureGenerator::align_result& r) {

    std::vector<std::string> output_rows;
    std::uint32_t num_rows = r.target_positions_pileup[0].size();
    output_rows.resize(num_rows);
    for (auto& row: output_rows) {
        row.reserve(r.target_positions_pileup.size());
    }
    
    for (std::uint32_t i = 0; i < r.target_positions_pileup.size(); i++) { // for each target position        
        auto& column = r.target_positions_pileup[i];
        for (std::uint32_t k = 0; k < num_rows; k++) { // move down the column
            output_rows[k].push_back(DECODER[column[k]]);                                              
        }
        
        auto& ins_columns = r.ins_positions_pileup[i];
        for (std::uint32_t j = 0; j < ins_columns.size() ; j++) { // for each ins column here
            auto& ins_column = ins_columns[j];                       
            for (std::uint32_t k = 0; k < num_rows; k++) { // move down the column
                output_rows[k].push_back(DECODER[ins_column[k]]);                   
            }
        }
           
    }
    //printing
    for (auto& row: output_rows) {
        std::cout << row << std::endl;
        
    }
    
    
}





       
} // namespace align_reads