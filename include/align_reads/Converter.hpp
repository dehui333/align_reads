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
 *
 */

namespace align_reads 
{
    // Package carrying feature matrices and auxiliary info
    // like which sequence is being aligned to & which sequences are aligning
    struct Data 
    {
        std::vector<PyObject*> Xs;
    };

    class AlignmentConverter
    {
        public:
            AlignmentConverter(MultiAlignment& alignment, std::uint32_t matrix_height, std::uint32_t matrix_width);

            // Insert input feature matrices 
            void produce_data();

        private:
            MultiAlignment* alignment_ptr;
            std::uint32_t matrix_width;
            std::uint32_t matrix_height;
            std::vector<std::pair<std::uint32_t, std::uint32_t>> width_idx_to_pos_idx;
            std::vector<std::vector<std::uint32_t>> segments_in_windows; 
            std::vector<std::pair<std::uint32_t, std::uint32_t>> start_of_windows;
            
            // may be inlined?
            void fill_row_from_alignment(PyObject* matrix, std::uint32_t window, 
                std::uint32_t row, std::uint32_t alignment_idx);

            // may be inlined?
            void fill_row_from_target(PyObject* matrix, std::uint32_t window, 
                std::uint32_t row);

            // Choose which segments to put into rows of each matrix, where they are 0 indexed.
            // Index 'number of alignments' will refer to fill the row with the target. 
            std::vector<std::vector<std::uint32_t>> choose_segments(std::uint32_t num_reserved_for_target, bool sample_target);



            
    };

} // namespace align_reads

#endif // ALIGN_READS_CONVERTER_HPP_
