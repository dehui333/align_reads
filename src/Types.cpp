#include <assert.h>

#include "edlib.h"

#include "align_reads/Types.hpp"

#define GAP_CHAR '_'

namespace align_reads
{
    alignment_segment::alignment_segment(std::string &query, std::string &target,
                                         std::uint32_t q_start, std::uint32_t t_start,
                                         EdlibAlignResult &result)
    {
        std::uint32_t aligned_len_on_target = result.endLocations[0] + 1 - result.startLocations[0];
        this->start_on_target = t_start + result.startLocations[0];
        this->end_on_target = t_start + result.endLocations[0] + 1;
        this->aligned_chars.reserve(aligned_len_on_target);
        this->ins_segments.resize(aligned_len_on_target + 1);
        std::uint32_t dest_index = 0;
        const char *current_q_char = query.c_str() + q_start;
        const char *current_t_char = target.c_str() + start_on_target;
        std::string ins_segment_buffer;
        ins_segment_buffer.reserve(1000);

        for (int i = 0; i < result.alignmentLength; i++)
        {
            switch (result.alignment[i])
            {
            case 0: // match
            {
                if (!ins_segment_buffer.empty())
                {
                    ins_segments[dest_index] = ins_segment_buffer;
                    ins_segment_buffer.clear();
                }
                dest_index++;
                this->aligned_chars.push_back(*current_q_char);
                current_q_char++;
                current_t_char++;
                break;
            }
            case 1: // ins
            {
                ins_segment_buffer.push_back(*(current_q_char++));
                break;
            }
            case 2: // del
            {
                assert(ins_segment_buffer.empty());
                dest_index++;
                this->aligned_chars.push_back(GAP_CHAR);
                current_t_char++;
                break;
            }
            case 3: // mismatch
            {
                if (!ins_segment_buffer.empty())
                {
                    ins_segments[dest_index] = ins_segment_buffer;
                    ins_segment_buffer.clear();
                }
                dest_index++;
                this->aligned_chars.push_back(*current_q_char);
                current_q_char++;
                current_t_char++;
                break;
            }
            }
        }
        if (!ins_segment_buffer.empty())
        {
            ins_segments[dest_index] = ins_segment_buffer;
        }
        edlibFreeAlignResult(result);
    }

    std::string &alignment_segment::get_aligned_chars()
    {
        return aligned_chars;
    }

    std::string &alignment_segment::get_ins_segment_at(int index)
    {
        return ins_segments[index + 1];
    }

}