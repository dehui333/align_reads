#ifndef ALIGN_READS_CONVERTER_HPP_
#define ALIGN_READS_CONVERTER_HPP_

#include <Python.h>
#include <utility>

#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#define NO_IMPORT_ARRAY
#include "numpy/arrayobject.h"

#include "align_reads/Aligner.hpp"


/*
 * Takes alignment objects and produce feature matrices.
 * - exposes inner implementation of some stuff in Aligner.
 */

namespace align_reads 
{
    // Package carrying feature matrices and auxiliary info
    // like which sequence is being aligned to & which sequences are aligning
    struct Data 
    {
        std::vector<PyObject*> Xs;
    };

    // Stores Additional info about the alignment/sequences 
    struct Info 
    {
        std::uint8_t hap_of_target;
        std::vector<std::uint8_t> hap_of_aligned;
    };

    // An 'adaptor' over a MultiAlginment object to produce feature matrices from it
    class AlignmentConverter
    {
        public:
            AlignmentConverter(MultiAlignment& alignment, std::uint32_t matrix_height, std::uint32_t matrix_width);
            AlignmentConverter(MultiAlignment& alignment, std::uint32_t matrix_height, std::uint32_t matrix_width, Info& info);

            // input feature matrices 
            Data produce_data(std::shared_ptr<thread_pool::ThreadPool> &pool, bool with_labels=false);

            



        private:
            MultiAlignment* alignment_ptr;
            Info* info_ptr;
            std::uint32_t matrix_width;
            std::uint32_t matrix_height;
            // a mapping from 'width index' to 'pos index'
            // width index: the column indices along the width of the alignment
            // pos index: An index which takes into account which position is aligned
            // to an existing base on the target and which are insertions etc
            std::vector<std::pair<std::uint32_t, std::uint32_t>> width_idx_to_pos_idx;

            // specifies which segments falls into which of the windows
            std::vector<std::vector<std::uint32_t>> segments_in_windows; 
            //std::vector<std::pair<std::uint32_t, std::uint32_t>> start_of_windows;
            
            // Fill a row of the matrix with bases from the aligned segment.
            // The 'window' parameter determines which part of the segment is taken.
            void fill_row_from_alignment(PyObject* matrix, std::uint32_t window, 
                std::uint32_t row, AlignmentSegment& segment);

            // Fill from the alignement target.
            void fill_row_from_target(PyObject* matrix, std::uint32_t window, 
                std::uint32_t row);

            // Choose which segments to put into rows of each matrix, where they are 0 indexed.
            // The output will specify the chosen segments to fill the matrices for each window.
            std::vector<std::vector<std::uint32_t>> choose_segments(std::uint32_t num_reserved_for_target, bool sample_target);

            // Produce matrices which consist of alignments of sequences aligned to the targert
            std::vector<PyObject*> produce_alignment_matrices(std::vector<std::vector<std::uint32_t>>& chosen, std::shared_ptr<thread_pool::ThreadPool> &pool);

            
    };

} // namespace align_reads

#endif // ALIGN_READS_CONVERTER_HPP_
