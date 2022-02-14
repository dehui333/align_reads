#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <utility>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "bioparser/paf_parser.hpp"
#include "feature_generator.hpp"
#include "types.hpp"

namespace align_reads {
    
FeatureGenerator::FeatureGenerator(const char* sequences_path, const char* overlaps_path) {
    std::string seq_path {sequences_path};
    std::string ov_path {overlaps_path};
    
    auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
    };
    
    if (is_suffix(seq_path, ".fasta")) {
        try { 
            auto p = bioparser::Parser<align_reads::Sequence>::Create<bioparser::FastaParser>(seq_path);
            auto s = p->Parse(-1);
            sequences.insert(sequences.end(), std::make_move_iterator(s.begin()), std::make_move_iterator(s.end()));
        } catch (const std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
        }            
              
    } else if (is_suffix(seq_path, ".fastq")) {
        try { 
            auto p = bioparser::Parser<align_reads::QSequence>::Create<bioparser::FastqParser>(seq_path);
            auto s = p->Parse(-1);
            sequences.insert(sequences.end(), std::make_move_iterator(s.begin()), std::make_move_iterator(s.end()));
        } catch (const std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
        }            

    } else {
        throw std::invalid_argument("[align_reads::FeatureGenerator::FeatureGenerator] Error: Invalid sequences file format.");
    }
    
    if (is_suffix(ov_path, ".paf")) {
        try { 
            auto p = bioparser::Parser<align_reads::Overlap>::Create<bioparser::PafParser>(ov_path);
            auto o = p->Parse(-1);
            overlaps.insert(overlaps.end(), std::make_move_iterator(o.begin()), std::make_move_iterator(o.end()));
        } catch (const std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
        }      
    } else {
        throw std::invalid_argument("[align_reads::FeatureGenerator::FeatureGenerator] Error: Invalid overlaps file format.");        
    }
}
       
}