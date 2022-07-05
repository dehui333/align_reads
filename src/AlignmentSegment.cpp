
#include "align_reads/AlignmentSegment.hpp"

#define GAP_CHAR '_'
#define PAD_CHAR '*' // Also defined in Converter.cpp

namespace align_reads
{
    //-----------------AlignmentSegment---------------------------

    AlignmentSegment::AlignmentSegment(std::string &query, std::uint32_t q_start,
                                       std::string &target, std::uint32_t t_start,
                                       EdlibAlignResult &result)
    {
        std::uint32_t aligned_len_on_target = result.endLocations[0] + 1 - result.startLocations[0];
        this->start_on_target = t_start + result.startLocations[0];
        this->end_on_target = t_start + result.endLocations[0];
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

    std::string &AlignmentSegment::get_aligned_chars()
    {
        return aligned_chars;
    }

    std::string &AlignmentSegment::get_ins_segment_at(int index)
    {
        return ins_segments[index + 1];
    }

    char AlignmentSegment::get_at_target_pos(std::uint32_t target_index, std::uint32_t ins_index)
    {
        // Nons ins positions
        if (ins_index == 0)
            return aligned_chars[target_index - start_on_target];

        // Ins positions
        if (ins_index > ins_segments[target_index - start_on_target + 1].size())
        {
            return GAP_CHAR;
        }

        return ins_segments[target_index - start_on_target + 1][ins_index - 1];
    }

    void AlignmentSegment::print(std::string &target)
    {
        // find out the 'width' - the total number of columns
        std::uint32_t width = aligned_chars.size();
        for (auto &s : ins_segments)
        {
            width += s.size();
        }

        // The alignment would be displayed in blocks
        std::uint16_t block_size = 100;
        std::uint32_t number_of_blocks = width % block_size == 0 ? width / block_size : width / block_size + 1;

        // One chunk from the target + one chunk from the aligned will form one block
        std::vector<std::string> chunks;
        chunks.resize(number_of_blocks * 2);
    }

}
