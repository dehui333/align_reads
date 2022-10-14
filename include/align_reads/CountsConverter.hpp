#ifndef COUNTS_CONVERTER_HPP_
#define COUNTS_CONVERTER_HPP_

#include <Python.h>
#include "align_reads/AlignCounter.hpp"

#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#define NO_IMPORT_ARRAY
#include "numpy/arrayobject.h"

namespace align_reads
{
    class CountsConverter
    {
    public:
        /* sequential production of matrices only.
        numpy array allocation is not thread safe.
        for parallelism, need to allocate the right number of matrices sequentially and determine
        beforehand the start and end counter (maybe using iterator) of each window  
        */
        static std::vector<PyObject*> get_counts_matrices(AlignCounter &counter, std::uint16_t window_length, std::uint32_t left_clip, std::uint32_t right_clip, std::uint32_t num_matrices, std::uint32_t alignment_length);


        //CountsConverter(std::uint16_t window_length) : window_length(window_length) {}

    private:
        //std::uint16_t window_length; 
        
    };

} // namespace align_reads

#endif // COUNTS_CONVERTER_HPP_