#include <cassert>
#include <iostream>

#include "align_reads/AlignmentSegment.hpp"

#define GAP_CHAR '_'
#define PAD_CHAR '*' // Also defined in Converter.cpp

namespace align_reads
{
    bool AlignmentSegment::AlignmentIterator::has_next()
    {
        return current_target_index <= segment->end_on_target;
    }

    aligned_pos AlignmentSegment::AlignmentIterator::next()
    {
        if (current_ins_index == 0)
        {
            aligned_pos r{segment->aligned_chars[current_target_index - segment->start_on_target], current_target_index, current_ins_index};
            if (current_ins_index == current_max_ins)
            {
                current_target_index++;
                current_ins_index = 0;
                current_max_ins = segment->ins_segments[current_target_index - segment->start_on_target + 1].size();
            }
            else
            {
                current_ins_index++;
            }
            return r;
        }
        else
        {
            aligned_pos r;
            if (current_ins_index > current_max_ins)
            {
                //should throw exception
                r.c = PAD_CHAR;
                r.target_index = current_target_index;
                r.ins_index = current_ins_index;
            }
            else
            {
                // aligned_pos r{segment->ins_segments[current_target_index - segment->start_on_target + 1][current_ins_index - 1], current_target_index, current_ins_index};
                r.c = segment->ins_segments[current_target_index - segment->start_on_target + 1][current_ins_index - 1];
                r.target_index = current_target_index;
                r.ins_index = current_ins_index;
            }

            if (current_ins_index >= current_max_ins)
            {
                current_target_index++;
                current_ins_index = 0;
                current_max_ins = segment->ins_segments[current_target_index - segment->start_on_target + 1].size();
            }
            else
            {
                current_ins_index++;
            }
            return r;
        }
    }

    AlignmentSegment::AlignmentIterator::AlignmentIterator(AlignmentSegment &segment, std::uint32_t start_t_idx, std::uint32_t start_i_idx)
        : segment(&segment), current_target_index(start_t_idx), current_ins_index(start_i_idx)
    {
        current_max_ins = segment.ins_segments[start_t_idx - segment.start_on_target + 1].size();
    }


    //-----------------AlignmentSegment---------------------------
    AlignmentSegment::AlignmentSegment(std::string &target)
    {
        this->aligned_chars = target;
        this->start_on_target = 0;
        this->end_on_target = target.size() - 1;
        this->ins_segments.resize(target.size() + 1);
    }

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

    AlignmentSegment::AlignmentIterator AlignmentSegment::iterator(std::uint32_t start_t_idx, std::uint32_t start_i_idx)
    {
        return {*this, start_t_idx, start_i_idx};
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
        for (auto &s : chunks)
        {
            s.reserve(block_size);
        }

        // Fill up the blocks
        std::uint32_t width_idx = 0;
        auto iter = this->iterator(this->start_on_target, 0);
        std::vector<std::pair<std::uint32_t, std::uint32_t>> start_of_blocks; // Record positions
        start_of_blocks.reserve(number_of_blocks);
        std::uint32_t block_idx = 0;
        while (iter.has_next())
        {
            block_idx = width_idx / block_size;
            auto pos = iter.next();
            char target_char = pos.ins_index == 0 ? target[pos.target_index] : GAP_CHAR;
            chunks[block_idx * 2].push_back(target_char);
            chunks[block_idx * 2 + 1].push_back(pos.c);
            if (width_idx % block_size == 0)
                start_of_blocks.emplace_back(pos.target_index, pos.ins_index);
            width_idx++;
        }

        for (block_idx = 0; block_idx < number_of_blocks; block_idx++)
        {
            std::cout << "Block: " << block_idx << std::endl;
            auto start = start_of_blocks[block_idx];
            std::cout << "Start index: " << start.first << ", " << start.second << std::endl;
            std::cout << chunks[block_idx * 2] << std::endl;
            std::cout << chunks[block_idx * 2 + 1] << std::endl
                      << std::endl;
        }
    }

}
