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
#include "edlib.h"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

namespace align_reads {

class Aligner {

private:    
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
    std::shared_ptr<thread_pool::ThreadPool> pool;
    ram::MinimizerEngine minimizer_engine;
    std::vector<std::uint32_t> id_to_pos_index;
   
       
    struct align_result {
        std::vector<std::vector<char>> target_columns;
        // which target position? -> which ins at that position? -> which row?
        std::vector<std::vector<std::vector<char>>> ins_columns;
        std::uint32_t width; 
        
        
        align_result() = default;
        align_result(align_result&& r) : target_columns(std::move(r.target_columns)), ins_columns(std::move(r.ins_columns)), 
            width(r.width) {};
        
        align_result& operator= (align_result&& r) {
            if (this != &r) {
                this->target_columns = std::move(r.target_columns);
                this->ins_columns = std::move(r.ins_columns);
                this->width = r.width;
            }
            return *this;    
        };        
        void print();
    };
    
           
    
    // align queries to target
    static align_result align_to_target_clip(std::vector<std::string>& queries, std::string& target,
        std::vector<std::pair<std::uint32_t, std::uint32_t>>& clips, std::vector<std::vector<std::uint32_t>>& inserters,
        std::vector<std::uint32_t>& ins_at_least2, bool has_hap, std::vector<EdlibAlignResult>& edlib_results);
    
    static align_result align_to_target_no_clip(std::vector<std::string>& queries, std::string& target, bool has_hap);
    
    // align queries to target, and also try to align the ins segments
    static align_result multi_align(std::vector<std::string>& queries, std::string& target,
        std::vector<std::pair<std::uint32_t, std::uint32_t>>& clips, std::vector<EdlibAlignResult>& edlib_results, bool has_hap=false);
           
    

public:
   
    
    Aligner(const char** sequences_file_paths, std::shared_ptr<thread_pool::ThreadPool>& pool,
        std::uint8_t kmer_len, std::uint8_t window_len, double freq);
     
    void find_true_overlaps();

    void find_RAM_overlaps();

    void within_each();

    void true_positives();

    void false_positives();

    void false_negatives();

    void run();
    
};



} // namespace align_reads
#endif  // ALIGN_READS_FEATURE_GENERATOR_HPP_
