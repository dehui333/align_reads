#ifndef ALIGN_READS_FEATURE_GENERATOR_HPP_
#define ALIGN_READS_FEATURE_GENERATOR_HPP_

#include <Python.h>
//extern "C" {
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#include "numpy/arrayobject.h"
//}

#include <atomic>
#include <memory> // unique_ptr
#include <utility>
#include <vector> 

#include "biosoup/nucleic_acid.hpp"
#include "ram/minimizer_engine.hpp"

namespace align_reads {

struct Data {
    std::vector<PyObject*> X;
    std::vector<PyObject*> Y;        
};

class Aligner {

private:    
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    ram::MinimizerEngine minimizer_engine;
    std::vector<std::uint32_t> id_to_pos_index;
    
    std::uint32_t start_of_other_phase = -1;
    
    std::vector<std::vector<std::unique_ptr<biosoup::NucleicAcid>>> haplotypes_sequences;
    std::vector<ram::MinimizerEngine> haplotypes_minimizer_engines;
    
    std::atomic<std::uint32_t> num_processed{0};
    
    
    struct seq_info {
        std::uint32_t id;
        bool reverse_complement; // is actually seq the reverse complement of what's in the alignment?
        seq_info(std::uint32_t id, bool reverse_comp) : id(id), reverse_complement(reverse_comp) {};
    };
    
    struct align_boundary {
		// refer to target indices
        std::uint32_t align_start; 
        std::uint32_t align_end;   // inclusive 
        bool contained;
        align_boundary(std::uint32_t start, std::uint32_t end, bool contained) : align_start(start), align_end(end), contained(contained) {};
    };
    
    struct align_result {
        std::vector<std::vector<char>> target_columns;
        // which target position? -> which ins at that position? -> which row?
        std::vector<std::vector<std::vector<char>>> ins_columns;
        std::uint32_t width; 
        std::vector<std::vector<std::uint32_t>> inserters; // for each position, contains positional index of queries that has ins there
        std::vector<std::uint32_t> ins_at_least2; // target positions where there are at least 2 reads with ins
        std::vector<align_boundary> align_boundaries; // where do the others align to the top sequence (target)
        
        
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
        
        void print();
    };
    
    struct align_overlapping_result {
        
        align_result alignment;
        std::vector<seq_info> infos;
        std::uint32_t target_id;
        
        align_overlapping_result() = default;
        align_overlapping_result(align_result&& alignment, std::vector<seq_info>&& infos, std::uint32_t& target_id) 
            : alignment(std::move(alignment)), infos(std::move(infos)), target_id(target_id) {};
        align_overlapping_result(align_overlapping_result&& r) : alignment(std::move(r.alignment)), infos(std::move(r.infos)),
            target_id(r.target_id) {};
        align_overlapping_result& operator= (align_overlapping_result&& r) {
            this->alignment = std::move(r.alignment);
            this->infos = std::move(r.infos);
            this->target_id = r.target_id;
            return *this;
        };
        
        Data produce_data(bool produce_labels=false, std::uint32_t start_of_other_phase=0);
    };
    
    
    
    align_overlapping_result align_overlapping_plus_haplotypes(std::unique_ptr<biosoup::NucleicAcid>& seq);
    
    // align queries to target
    static align_result align_to_target(std::vector<std::string>& queries, std::string& target,
        bool clip_query, std::vector<std::pair<std::uint32_t, std::uint32_t>>* pads=nullptr);
    
    align_overlapping_result align_overlapping(std::unique_ptr<biosoup::NucleicAcid>& seq);
    
    // align queries to target, and also try to align the ins segments
    static align_result pseudoMSA(std::vector<std::string>& queries, std::string& target,
        std::vector<std::pair<std::uint32_t, std::uint32_t>>& pads);
        


public:
   
    
    Aligner(const char** sequences_file_paths, std::uint32_t num_threads,
        std::uint8_t kmer_len, std::uint8_t window_len, double freq, const char** haplotypes_path=nullptr);
     
    Data next();
    
    bool has_next();
    
};



} // namespace align_reads
#endif  // ALIGN_READS_FEATURE_GENERATOR_HPP_
