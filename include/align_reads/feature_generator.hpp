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
    
    std::uint32_t start_of_other_phase = -1;
    
    std::vector<std::vector<std::unique_ptr<biosoup::NucleicAcid>>> haplotypes_sequences;
    std::vector<ram::MinimizerEngine> haplotypes_minimizer_engines;
    std::vector<std::vector<std::uint32_t>> haplotypes_id_to_pos_index;
    
    struct read_info {
        std::uint32_t id;
        bool reverse_complement;
        read_info(std::uint32_t id, bool reverse_comp) : id(id), reverse_complement(reverse_comp) {};
    };
    
    struct align_boundary {
        std::uint32_t align_start; 
        std::uint32_t align_end;   // inclusive 
        align_boundary(std::uint32_t start, std::uint32_t end) : align_start(start), align_end(end) {};
    };
    
    struct align_result {
        std::vector<std::vector<char>> target_columns;
        // which target position? -> which ins at that position? -> which row?
        std::vector<std::vector<std::vector<char>>> ins_columns;
        std::uint32_t width;
        std::vector<std::vector<std::uint32_t>> inserters; // contains positional index of queries
        std::vector<std::uint32_t> ins_at_least2;
        std::vector<align_boundary> align_boundaries;
        
        align_result() = default;
        align_result(align_result&& r) : target_columns(std::move(r.target_columns)), ins_columns(std::move(r.ins_columns)), 
            width(r.width), inserters(std::move(r.inserters)), ins_at_least2(std::move(r.ins_at_least2)),
            align_boundaries(std::move(r.align_boundaries)) {};
        
        align_result& operator= (align_result&& r) {
            if (this != &r) {
                this->target_columns = std::move(r.target_columns);
                this->ins_columns = std::move(r.ins_columns);
                this->inserters = std::move(r.inserters);
                this->ins_at_least2 = std::move(r.ins_at_least2);
                this->width = r.width;
                this->align_boundaries = std::move(r.align_boundaries);
            }
            return *this;    
        };
    };
    
    struct align_overlapping_result {
        align_result alignment;
        std::vector<read_info> infos;
        std::uint32_t target_id;
        
        align_overlapping_result() = default;
        align_overlapping_result(align_result&& alignment, std::vector<read_info>&& infos, std::uint32_t& target_id) 
            : alignment(std::move(alignment)), infos(std::move(infos)), target_id(target_id) {};
        align_overlapping_result(align_overlapping_result&& r) : alignment(std::move(r.alignment)), infos(std::move(r.infos)),
            target_id(r.target_id) {};
        align_overlapping_result& operator= (align_overlapping_result&& r) {
            this->alignment = std::move(r.alignment);
            this->infos = std::move(r.infos);
            this->target_id = r.target_id;
            return *this;
        };
    };
        
    
    align_overlapping_result align_overlapping(std::unique_ptr<biosoup::NucleicAcid>& seq);
    
    // align queries to target
    static align_result align_to_target(std::vector<std::string>& queries, std::string& target, bool clip_query);
    
    // align queries to target, and also try to align the ins segments
    static align_result pseudoMSA(std::vector<std::string>& queries, std::string& target);
    
    static void print_align(align_result& r);
public:
    
    FeatureGenerator(const char** sequences_file_paths, std::uint32_t num_threads,
        std::uint8_t kmer_len, std::uint8_t window_len, double freq, const char** haplotypes_path=nullptr);
    
    
};



} // namespace align_reads
#endif  // ALIGN_READS_FEATURE_GENERATOR_HPP_
