#ifndef ALIGN_COUNTER_HPP_
#define ALIGN_COUNTER_HPP_

#include <string>
#include <vector>

#include "edlib.h"

#include "align_reads/Aligner.hpp"

namespace align_reads
{

    constexpr static std::uint8_t ENCODER[] = {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 5, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 0, 255, 1, 255, 255,
        255, 2, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 3, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 4, 255, 0, 255, 1,
        255, 255, 255, 2, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 3, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255};
    constexpr static char DECODER[] = {
        'A', 'C', 'G', 'T', '_', '*'};

    class AlignCounter
    {
    public:
        // frees the EdlibAlignResult
        AlignCounter(std::string &target_string,
                     std::vector<clipped_alignment<EdlibAlignResult>> &alignments);
        inline std::uint32_t max_ins_at(std::uint32_t target_idx)
        {
            return counts[target_idx].size() - 1;
        }

        inline std::uint16_t base_count_at(std::uint32_t target_idx, std::uint16_t ins_idx, std::uint8_t base_idx)
        {
            return counts[target_idx][ins_idx][base_idx];
        }

    private:
        std::vector<std::vector<std::vector<std::uint16_t>>> counts;
    };

} // namespace align_reads

#endif // ALIGN_COUNTER_HPP_