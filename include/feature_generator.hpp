#ifndef ALIGN_READS_FEATURE_GENERATOR_HPP_
#define ALIGN_READS_FEATURE_GENERATOR_HPP_

#include <memory> // unique_ptr
#include <vector> 

#include "types.hpp"

namespace align_reads {

class FeatureGenerator {

private:    
    std::vector<std::unique_ptr<align_reads::Sequence>> sequences;
    std::vector<std::unique_ptr<align_reads::Overlap>> overlaps;
    
public:
    
    FeatureGenerator(const char* sequences_file_path, const char* overlaps_file_path);
      
};



} // namespace align_reads
#endif  // ALIGN_READS_FEATURE_GENERATOR_HPP_
