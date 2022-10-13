#include <iostream>

#include "align_reads/AlignCounter.hpp"

#define A_INDEX 0
#define C_INDEX 1
#define G_INDEX 2
#define T_INDEX 3
#define DEL_INDEX 4

namespace align_reads
{

    inline void initialize_slot(std::vector<std::vector<std::uint16_t>> &counts_at_pos, std::uint32_t slot_index, std::uint32_t& alignment_length)
    {
        if (counts_at_pos.size() < slot_index + 1)
        {
            alignment_length++;
            counts_at_pos.resize(slot_index + 1);
            counts_at_pos[slot_index].resize(5, 0);
        }
    }

    inline void update_counts_and_stats(std::vector<std::vector<std::vector<std::uint16_t>>> &counts,
                                        std::vector<std::vector<std::vector<float>>>& stats,
                                        clipped_alignment<EdlibAlignResult> &alignment, std::uint32_t& alignment_length)

    {
        auto &result = alignment.result;
        const char *current_q_char = alignment.clipped_query.c_str();
        int start_t_idx = static_cast<int>(alignment.t_start) + result.startLocations[0] - 1;
        int t_idx = start_t_idx;
        std::uint32_t ins_idx = 1;
        for (int i = 0; i < result.alignmentLength; i++)
        {
            switch (result.alignment[i])
            {
            case 0: // match
            case 3: // mismatch
            {

                counts[++t_idx][0][ENCODER[static_cast<std::uint8_t>(*current_q_char)]]++;

                current_q_char++;
                ins_idx = 1;
                if (static_cast<std::uint32_t>(t_idx) == counts.size() - 1)
                {
                    // ignore ins at the end
                    i = result.alignmentLength;
                }
                break;
            }
            case 1: // ins
            {
                if (t_idx == start_t_idx)
                {
                    // ignore ins at start
                    current_q_char++;
                    break;
                }
                std::vector<std::vector<std::uint16_t>> &counts_at_pos = counts[t_idx];
                initialize_slot(counts_at_pos, ins_idx, alignment_length);
                counts_at_pos[ins_idx++][ENCODER[static_cast<std::uint8_t>(*current_q_char)]]++;
                current_q_char++;
                break;
            }
            case 2: // del
            {
                counts[++t_idx][0][4]++;
                break;
            }
            }
        }
        edlibFreeAlignResult(result);
    }

    AlignCounter::AlignCounter(std::string &target_string,
                               std::vector<clipped_alignment<EdlibAlignResult>> &alignments)
    {
        // initialize and add counts from target
        counts.resize(target_string.size());
        //stats.resize(target_string.size());
        alignment_length = target_string.size();
        std::uint32_t i;
        for (i = 0; i < target_string.size(); i++)
        {
            auto &counts_at_pos = counts[i];
            counts_at_pos.resize(1);
            counts_at_pos[0].resize(5, 0);
            counts_at_pos[0][ENCODER[static_cast<std::uint8_t>(target_string[i])]]++;

            /*
            auto &counts_at_pos = counts[i];
            counts_at_pos.resize(1);
            counts_at_pos[0].resize(5, 0);
            counts_at_pos[0][ENCODER[static_cast<std::uint8_t>(target_string[i])]] = 1;
            */
        }

        // Add in counts from aligned queries
        for (auto &alignment : alignments)
        {
            update_counts_and_stats(counts, stats, alignment, alignment_length);
        }

        // add in del for those without bases at some ins pos
        for (auto &counts_at_pos : counts)
        {
            auto &counter0 = counts_at_pos[0];
            std::uint16_t num_support = counter0[0] + counter0[1] + counter0[2] + counter0[3] + counter0[4];
            for (i = 1; i < counts_at_pos.size(); i++)
            {
                auto &counter = counts_at_pos[i];
                counter[4] += (num_support - counter[0] - counter[1] - counter[2] - counter[3]);
            }
        }
    }

} // namespace align_reads