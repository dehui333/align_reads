#ifndef ALIGN_READS_FEATURE_GENERATOR_HPP_
#define ALIGN_READS_FEATURE_GENERATOR_HPP_

#include <memory> // unique_ptr
#include <vector> 

#include "biosoup/nucleic_acid.hpp"
#include "ram/minimizer_engine.hpp"

namespace align_reads {

class FeatureGenerator {

private:    
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    ram::MinimizerEngine minimizer_engine;
    
public:
    
    FeatureGenerator(const char* sequences_file_path, std::uint32_t num_threads,
        std::uint8_t kmer_len, std::uint8_t window_len, double freq);
    
    //std::vector<std::uint32_t> find_overlapping(std::unique_tr<biosoup::NucleicAcid> seq);   
};



} // namespace align_reads
#endif  // ALIGN_READS_FEATURE_GENERATOR_HPP_
