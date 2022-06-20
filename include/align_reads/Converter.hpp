#ifndef ALIGN_READS_CONVERTER_HPP_
#define ALIGN_READS_CONVERTER_HPP_

#include "align_reads/Aligner.hpp"

/*
 * Takes alignment objects and produce feature matrices.
 *
 */

namespace align_reads 
{
    class AlignmentConverter
    {
        public:
            AlignmentConverter(MultiAlignment& alignment, std::uint32_t matrix_width, std::uint32_t matrix_height);
        private:
            MultiAlignment* alignment_ptr;
            std::uint32_t matrix_width;
            std::uint32_t matrix_height;
            std::vector<std::vector<std::uint32_t>> segments_in_windows; 
            
    };

} // namespace align_reads

#endif // ALIGN_READS_CONVERTER_HPP_
