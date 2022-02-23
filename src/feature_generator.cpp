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

FeatureGenerator::align_result FeatureGenerator::align(std::unique_ptr<biosoup::NucleicAcid>& target) {
    std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, true); // find all overlaps with target 
    auto sort_by_id = [] (const biosoup::Overlap& o1, const biosoup::Overlap& o2) {
        return o1.rhs_id < o2.rhs_id;        
    };
    std::vector<biosoup::Overlap*> unique_overlaps; // sort overlaps by id 
    std::sort(overlaps.begin(), overlaps.end(), sort_by_id);
    std::uint64_t last_id = -1;
    // remove duplicates
    for (auto& o : overlaps) {
        if (o.rhs_id != last_id) {
            unique_overlaps.push_back(&o);
            last_id = o.rhs_id;
        }
    }

    std::string target_underlying_string = target->InflateData();          
    const char* target_underlying = target_underlying_string.c_str();
    //start to fill in align results
    align_result result;
    result.infos.resize(unique_overlaps.size() + 1);
    result.target_positions_pileup.resize(target->inflated_len);
    result.ins_positions_pileup.resize(target->inflated_len);
    for (std::uint32_t i = 0; i < result.target_positions_pileup.size(); i++) {
        auto& v = result.target_positions_pileup[i];
        v.resize(unique_overlaps.size() + 1, ENCODER['_']);
        v[0] = ENCODER[static_cast<std::uint8_t>(target_underlying[i])];
    }
    result.infos[0].id = target->id;
    result.infos[0].same_strand = true;
    result.infos[0].align_start = 0;
    result.infos[0].align_end = target->inflated_len;
    result.width = target->inflated_len; // will be incremented
    for (auto& o : unique_overlaps) {
        std::uint32_t pos_index = id_to_pos_index[o->rhs_id];
        auto& query_seq = sequences[pos_index];
        if (!o->strand) {
            query_seq->ReverseAndComplement();
        }
        EdlibAlignResult align = edlibAlign(query_seq->InflateData().c_str(), query_seq->inflated_len,
            target_underlying, target->inflated_len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        std::uint32_t query_index = 0; // next query position
        std::uint32_t target_index = 0; // next target position
        std::uint32_t ins_index = 0; // next ins index for the current segment
        
        for (int i = 0; i < align.alignmentLength; i++) {
            switch (align.alignment[i]) {
                case 0:  // match
                    break;
                case 1:  // ins
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
    
    
    
    return result;    
}

void FeatureGenerator::print_align(FeatureGenerator::align_result& r) {
    /*
     struct read_info {
        std::uint32_t id;
        bool same_strand;
        std::uint32_t align_start;
        std::uint32_t align_end;
    };
    
    struct align_result {
        std::vector<read_info> infos;
        std::vector<std::vector<std::uint8_t>> target_positions_pileup;
        // which target position? -> which ins at that position? -> which row?
        std::vector<std::vector<std::vector<std::uint8_t>>> ins_positions_pileup;
    };
    */
    std::vector<std::string> output_rows;
    output_rows.resize(r.infos.size());
    for (auto& row: output_rows) {
        row.reserve(r.target_positions_pileup.size());
    }
    
    for (std::uint32_t i = 0; i < r.target_positions_pileup.size(); i++) { // for each target position        
        auto& column = r.target_positions_pileup[i];
        for (std::uint32_t k = 0; k < r.infos.size(); k++) { // move down the column
            output_rows[k].push_back(DECODER[column[k]]);                                              
        }
        
        auto& ins_columns = r.ins_positions_pileup[i];
        for (std::uint32_t j = 0; j < ins_columns.size() ; j++) { // for each ins column here
            auto& ins_column = ins_columns[j];                       
            for (std::uint32_t k = 0; k < r.infos.size(); k++) { // move down the column
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