#ifndef ALIGN_READS_FEATURE_GENERATOR_HPP_
#define ALIGN_READS_FEATURE_GENERATOR_HPP_

#include <memory> // unique_ptr
#include <utility>
#include <vector> 

#include "biosoup/nucleic_acid.hpp"
#include "ram/minimizer_engine.hpp"

namespace align_reads {

class FeatureGenerator {

private:    
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    ram::MinimizerEngine minimizer_engine;
    std::vector<std::uint32_t> id_to_pos_index;
    
    struct read_info {
        std::uint32_t id;
        bool same_strand;
        std::uint32_t align_start;
        std::uint32_t align_end;
    };
    
    struct align_result {
        std::vector<std::vector<char>> target_positions_pileup;
        // which target position? -> which ins at that position? -> which row?
        std::vector<std::vector<std::vector<char>>> ins_positions_pileup;
        std::uint32_t width;
        std::vector<std::vector<std::uint32_t>> inserters; // contains positional indices of queries
    };
    
    struct align_overlapping_result {
        align_result alignment;
        read_info info;
    };
        
    
    align_overlapping_result align_overlapping(std::unique_ptr<biosoup::NucleicAcid>& seq);
    
    // align queries to target
    static align_result align_to_target(std::vector<std::string>& queries, std::string& target, bool clip_query);
    
    // align queries to target, and also try to align the ins segments
    static align_result pseudoMSA(std::vector<std::string>& queries, std::string& target);
    
    static void print_align(align_result& r);
public:
    
    FeatureGenerator(const char* sequences_file_path, std::uint32_t num_threads,
        std::uint8_t kmer_len, std::uint8_t window_len, double freq);
    
    
};



} // namespace align_reads
#endif  // ALIGN_READS_FEATURE_GENERATOR_HPP_
