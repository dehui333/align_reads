#ifndef ALIGN_READS_FEATURE_GENERATOR_HPP_
#define ALIGN_READS_FEATURE_GENERATOR_HPP_

#include <memory> // unique_ptr
#include <unordered_set>
#include <vector> 

#include "biosoup/nucleic_acid.hpp"
#include "ram/minimizer_engine.hpp"

namespace align_reads {

class FeatureGenerator {

private:    
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    ram::MinimizerEngine minimizer_engine;
    std::vector<std::uint32_t> id_to_pos_index;
    
    struct reads_distribution {
        std::vector<std::unordered_set<std::uint32_t>> reads_by_segment;
        std::unordered_set<std::uint32_t> all_reads;
    };    
    
    reads_distribution distribute_reads(std::unique_ptr<biosoup::NucleicAcid>& seq);   
    void align(std::unique_ptr<biosoup::NucleicAcid>& seq);
public:
    
    FeatureGenerator(const char* sequences_file_path, std::uint32_t num_threads,
        std::uint8_t kmer_len, std::uint8_t window_len, double freq);
    
    
};



} // namespace align_reads
#endif  // ALIGN_READS_FEATURE_GENERATOR_HPP_
