#include <iostream>
#include <utility>

#include "align_reads/AlignCounter.hpp"

#define A_INDEX 0
#define C_INDEX 1
#define G_INDEX 2
#define T_INDEX 3
#define GAP_INDEX 4

namespace align_reads
{

    inline void initialize_slot(std::vector<std::vector<std::uint16_t>> &counts_at_pos, std::vector<std::vector<float>> &stats_at_pos, std::uint32_t slot_index, std::uint32_t &alignment_length)
    {
        if (counts_at_pos.size() < slot_index + 1)
        {
            if (slot_index > 0)
            {
                alignment_length++;
            }
            
            counts_at_pos.resize(slot_index + 1);
            counts_at_pos[slot_index].resize(5, 0);

            stats_at_pos.resize(slot_index + 1);
            stats_at_pos[slot_index].resize(5, 0);
        }
    }

    inline void update_counts_and_stats_one(
        std::vector<std::vector<std::uint16_t>> &counts_at_pos,
        std::vector<std::vector<float>> &stats_at_pos,
        std::uint16_t ins_idx,
        std::uint8_t base_idx,
        float identity_score)
    { 
        //counts_at_pos[ins_idx][base_idx]++;
        //stats_at_pos[ins_idx][base_idx] +=

        std::uint16_t& count = counts_at_pos[ins_idx][base_idx];
        float &stat = stats_at_pos[ins_idx][base_idx];
        count++;
        stat += (identity_score - stat)/count; // maintain average 
    }

    inline void update_counts_and_stats(std::vector<std::vector<std::vector<std::uint16_t>>> &counts,
                                        std::vector<std::vector<std::vector<float>>> &stats,
                                        clipped_alignment<EdlibAlignResult> &alignment, std::uint32_t &alignment_length,
                                        std::vector<std::vector<std::pair<std::uint16_t, float>>> &aux)

    {
        // Assumes infix alignment to avoid troublesome cases

        auto &result = alignment.result;
        const char *current_q_char = alignment.clipped_query.c_str();
        int start_t_idx = static_cast<int>(alignment.t_start) + result.startLocations[0] - 1;
        int t_idx = start_t_idx;
        std::uint32_t ins_idx = 1;
        std::uint32_t num_aligned_at_pos = 0;

        float identity_score = 1 - (float)result.editDistance / alignment.clipped_query.size();
        bool first_match_target = true;

        for (int i = 0; i < result.alignmentLength; i++)
        {
            switch (result.alignment[i])
            {
            case 0: // match
            case 3: // mismatch
            {
                if (first_match_target)
                {
                    first_match_target = false;
                }
                else
                {
                    aux[t_idx].emplace_back(num_aligned_at_pos, identity_score);
                }
                t_idx++;
                update_counts_and_stats_one(counts[t_idx], stats[t_idx], 0, ENCODER[static_cast<std::uint8_t>(*current_q_char)], identity_score);
                //counts[++t_idx][0][ENCODER[static_cast<std::uint8_t>(*current_q_char)]]++;
                current_q_char++;
                ins_idx = 1;
                num_aligned_at_pos = 1;
                if (static_cast<std::uint32_t>(t_idx) == counts.size() - 1)
                {
                    // ignore ins to the end of target
                    i = result.alignmentLength;
                }
                break;
            }
            case 1: // ins
            {
                if (t_idx == start_t_idx)
                {
                    // ignore ins to the start of target
                    current_q_char++;
                    break;
                }
                num_aligned_at_pos++;
                std::vector<std::vector<std::uint16_t>> &counts_at_pos = counts[t_idx];
                std::vector<std::vector<float>> &stats_at_pos = stats[t_idx];
                initialize_slot(counts_at_pos, stats_at_pos, ins_idx, alignment_length);
                update_counts_and_stats_one(counts_at_pos, stats_at_pos, ins_idx++, ENCODER[static_cast<std::uint8_t>(*current_q_char)], identity_score);
                //counts_at_pos[ins_idx++][ENCODER[static_cast<std::uint8_t>(*current_q_char)]]++;
                current_q_char++;
                break;
            }
            case 2: // del
            {
                // settle last chunk
                aux[t_idx].emplace_back(num_aligned_at_pos, identity_score);
                t_idx++;
                num_aligned_at_pos = 0;
                // counts[++t_idx][0][4]++; taken care of by aux
                break;
            }
            }
        }
        aux[t_idx].emplace_back(num_aligned_at_pos, identity_score);
        edlibFreeAlignResult(result);
    }

    AlignCounter::AlignCounter(std::string &target_string,
                               std::vector<clipped_alignment<EdlibAlignResult>> &alignments)
    {
        // initialize and add counts from target
        counts.resize(target_string.size());
        stats.resize(target_string.size());
        alignment_length = target_string.size();
        std::uint32_t i;
        for (i = 0; i < target_string.size(); i++)
        {
            auto &counts_at_pos = counts[i];
            //counts_at_pos.resize(1);
            //counts_at_pos[0].resize(5, 0);
            //counts_at_pos[0][ENCODER[static_cast<std::uint8_t>(target_string[i])]]++;

            auto &stats_at_pos = stats[i];
            //stats_at_pos.resize(1);
            //stats_at_pos[0].resize(5, 0);
            //stats_at_pos[0][ENCODER[static_cast<std::uint8_t>(target_string[i])]] += 1;
            initialize_slot(counts_at_pos, stats_at_pos, 0, alignment_length);
            update_counts_and_stats_one(counts_at_pos, stats_at_pos, 0, ENCODER[static_cast<std::uint8_t>(target_string[i])], 1);
            /*
            auto &counts_at_pos = counts[i];
            counts_at_pos.resize(1);
            counts_at_pos[0].resize(5, 0);
            counts_at_pos[0][ENCODER[static_cast<std::uint8_t>(target_string[i])]] = 1;
            */
        }

        // For each base on target, store a vector of
        // the identity score of aligned sequence + the number of bases aligned there
        // -> for identity score of gap in ins
        std::vector<std::vector<std::pair<std::uint16_t, float>>> aux;
        aux.resize(target_string.size());

        // Add in counts from aligned queries
        for (auto &alignment : alignments)
        {
            update_counts_and_stats(counts, stats, alignment, alignment_length, aux);
        }
        /*
        for (auto i = 0; i < aux.size(); i++)
        {
            std::cout << "t_idx " << i << std::endl;
            for (auto &p : aux[i])
            {
                std::cout << p.first << std::endl;
            }
        }*/

        // add in del for those without bases at some ins pos
        /*
        for (auto &counts_at_pos : counts)
        {
            auto &counter0 = counts_at_pos[0];
            std::uint16_t num_support = counter0[0] + counter0[1] + counter0[2] + counter0[3] + counter0[4];
            for (i = 1; i < counts_at_pos.size(); i++)
            {
                auto &counter = counts_at_pos[i];
                counter[4] += (num_support - counter[0] - counter[1] - counter[2] - counter[3]);
            }
        }*/

        i = 0;
        std::uint32_t j;
        for (auto &pairs_at_pos : aux)
        {
            auto &counts_at_pos = counts[i];
            auto &stats_at_pos = stats[i];
            for (auto &pair : pairs_at_pos)
            {
                // pair.first is the number of bases aligned to this target position by this
                // query
                for (j = pair.first; j < counts_at_pos.size(); j++)
                {
                    update_counts_and_stats_one(counts_at_pos, stats_at_pos, j, 4, pair.second);
                    //counts_at_pos[j][4]++;
                }
            }
            // add gap count for target
            for (j = 1; j < counts_at_pos.size(); j++)
            {
                //counts_at_pos[j][4]++;
                update_counts_and_stats_one(counts_at_pos, stats_at_pos, j, 4, 1);
            }
            i++;
        }
    }

} // namespace align_reads