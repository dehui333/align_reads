#ifndef ALIGN_READS_FEATURE_GENERATOR_HPP_
#define ALIGN_READS_FEATURE_GENERATOR_HPP_

#include <memory> // unique_ptr
#include <utility>
#include <vector> 

#include "biosoup/nucleic_acid.hpp"
#include "ram/minimizer_engine.hpp"

namespace align_reads {
    
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
        std::vector<read_info> infos;
        std::vector<std::vector<std::uint8_t>> target_positions_pileup;
        // which target position? -> which ins at that position? -> which row?
        std::vector<std::vector<std::vector<std::uint8_t>>> ins_positions_pileup;
        std::uint32_t width;
    };
        
        
        
    void print_align(align_result& r);
    
    align_result align(std::unique_ptr<biosoup::NucleicAcid>& seq);
public:
    
    FeatureGenerator(const char* sequences_file_path, std::uint32_t num_threads,
        std::uint8_t kmer_len, std::uint8_t window_len, double freq);
    
    
};



} // namespace align_reads
#endif  // ALIGN_READS_FEATURE_GENERATOR_HPP_
