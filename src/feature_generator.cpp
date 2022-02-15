#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>
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

// Find the id of reads which overlaps with seq
std::vector<std::uint32_t> FeatureGenerator::find_overlapping(std::unique_ptr<biosoup::NucleicAcid>& seq) {
    std::vector<std::uint32_t> overlapping_reads;
    std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(seq, true, false, true);
    overlapping_reads.reserve(overlaps.size());
    for (auto& o: overlaps) {
        overlapping_reads.push_back(o.rhs_id);
    }
    return overlapping_reads;
}

void FeatureGenerator::align(std::unique_ptr<biosoup::NucleicAcid>& seq) {
    std::vector<std::uint32_t> overlapping_reads = find_overlapping(seq);
    for (auto& id: overlapping_reads) {
        std::uint32_t pos_index = id_to_pos_index[id];
        EdlibAlignResult result = edlibAlign(sequences[pos_index]->InflateData().c_str(), sequences[pos_index]->inflated_len,
            seq->InflateData().c_str(), seq->inflated_len, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
        
        //do sth
        
        edlibFreeAlignResult(result); 
    }
    
        
}






       
}