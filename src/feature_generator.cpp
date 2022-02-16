#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <stdexcept>
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
constexpr static std::uint16_t MATRIX_COL_NUM = 256;
constexpr static std::uint16_t MATRIX_ROW_NUM = 256;

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

// Allocate overlapping reads to each segment of seq e.g. pos 1-> pos k, pos k + 1-> pos 2k ...
// Reads are allocated preferentially according to 'how contained' they are in seq
// After allocation, 
// If a segment has id count > MATRIX_ROW_NUM, its reads are filtered, preferring those that overlap more positions on seq
// If a segment has id count < MATRIX_ROW_NUM, do nothing for now, 
// wait till alignment to see if there are more actual overlaps, if still not enough, duplicate 
FeatureGenerator::reads_distribution FeatureGenerator::distribute_reads(std::unique_ptr<biosoup::NucleicAcid>& seq) {    
    std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(seq, true, false, true); // all overlaps with seq    
    std::uint16_t num_segments = seq->inflated_len / MATRIX_COL_NUM + (seq->inflated_len % MATRIX_COL_NUM != 0);
    std::vector<std::vector<biosoup::Overlap*>> overlaps_by_segment(num_segments);
    FeatureGenerator::reads_distribution distribution;
    
    // recalculate scores as the estimated number of non overlap positions on the query 
    std::uint32_t not_covered_left_query;
    std::uint32_t not_covered_right_query;
    std::uint32_t not_covered_left_target;
    std::uint32_t not_covered_right_target;
    int diff1;
    int diff2;
    for (auto& o: overlaps) {
        not_covered_left_query = o.lhs_begin;
        not_covered_right_query = MATRIX_COL_NUM * num_segments - o.lhs_end; // pad right with unknowns
        not_covered_left_target = o.rhs_begin;
        not_covered_right_target = sequences[id_to_pos_index[o.rhs_id]]->inflated_len - o.rhs_end;
        
        if (o.strand) {
            diff1 = not_covered_left_target - not_covered_left_query;
            diff2 = not_covered_right_target - not_covered_right_query;
            diff1 = diff1 >= 0 ? diff1 : 0;
            diff2 = diff2 >= 0 ? diff2 : 0;
        } else {
            diff1 = not_covered_left_target - not_covered_right_query;
            diff2 = not_covered_right_target - not_covered_left_query;
            diff1 = diff1 >= 0 ? diff1 : 0;
            diff2 = diff2 >= 0 ? diff2 : 0;
        }
        o.score = static_cast<std::uint32_t>(diff1 + diff2);
    }
    // sort
    auto least_not_covered_first = [] (const biosoup::Overlap& o1, const biosoup::Overlap& o2) {
        return o1.score < o2.score;
    };
    std::sort(overlaps.begin(), overlaps.end(), least_not_covered_first);
    
    // distribute 
    for (auto& o: overlaps) {
        std::uint32_t first_segment = o.lhs_begin / MATRIX_COL_NUM;
        std::uint32_t last_segment = (o.lhs_end - 1)/ MATRIX_COL_NUM;
        distribution.all_reads.insert(o.rhs_id);
        for (uint32_t i = first_segment; i <= last_segment; i++) {
            biosoup::Overlap* p = &o;
            overlaps_by_segment[i].push_back(p);
            distribution.reads_by_segment[i].insert(o.rhs_id);
        }
    }
    
    return distribution;
}

void FeatureGenerator::align(std::unique_ptr<biosoup::NucleicAcid>& seq) {
    FeatureGenerator::reads_distribution distribution = distribute_reads(seq);
    for (auto& id: distribution.all_reads) {
        std::uint32_t pos_index = id_to_pos_index[id];
        EdlibAlignResult result = edlibAlign(sequences[pos_index]->InflateData().c_str(), sequences[pos_index]->inflated_len,
            seq->InflateData().c_str(), seq->inflated_len, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        
        //do sth
        
        edlibFreeAlignResult(result); 
    }
    
        
}






       
}