#include <assert.h>

#include "align_reads/Converter.hpp"
#include "edlib.h"
#include "thread_pool/thread_pool.hpp"

#include "align_reads/Aligner.hpp"
#include "align_reads/Utilities.hpp"

namespace align_reads
{
    //------------------inlines----------------------
    /*inline std::uint32_t reduce_clip(std::uint32_t x)
    {
        if (x >= 20) return x - 20;
        return x;
    }*/
    //-----------------free---------------------------
    std::vector<EdlibAlignResult> get_edlib_results(std::vector<EdlibTask> &tasks,
                                                    std::shared_ptr<thread_pool::ThreadPool> &pool)
    {
        std::vector<std::future<EdlibAlignResult>> futures;
        std::vector<EdlibAlignResult> results;
        futures.reserve(tasks.size());
        results.reserve(tasks.size());

        for (auto &t : tasks)
        {
            futures.emplace_back(pool->Submit(
                edlibAlign, t.query_start,
                t.query_len, t.target_start,
                t.target_len, edlibNewAlignConfig(-1, t.mode, t.task, NULL, 0)));
        }

        for (auto &f : futures)
        {
            results.push_back(f.get());
        }
        return results;
    }

    EdlibAlignResult get_edlib_result(const char *q_start, std::uint32_t q_len,
                                      const char *t_start, std::uint32_t t_len, EdlibAlignMode mode, EdlibAlignTask task)
    {
        return edlibAlign(q_start, q_len, t_start, t_len, edlibNewAlignConfig(-1, mode, task, NULL, 0));
    }

    AlignmentSegment get_alignment_segment(std::string &query, std::uint32_t q_start, std::uint32_t q_len,
                                           std::string &target, std::uint32_t t_start, std::uint32_t t_len,
                                           EdlibAlignMode mode, EdlibAlignTask task)
    {
        auto result = edlibAlign(query.c_str() + q_start, q_len,
                                 target.c_str() + t_start, t_len, edlibNewAlignConfig(-1, mode, task, NULL, 0));
        return {query, q_start, target, t_start, result};
    }

    inline EdlibAlignResult get_edlib_result_(const char *q_start, std::uint32_t q_len,
                                              const char *t_start, std::uint32_t t_len, EdlibAlignMode mode, EdlibAlignTask task)
    {
        return edlibAlign(q_start, q_len, t_start, t_len, edlibNewAlignConfig(-1, mode, task, NULL, 0));
    }
    inline float iden_from_align_result(EdlibAlignResult &result)
    {
        std::uint32_t num_match = 0;
        std::uint32_t num_mismatch = 0;

        for (int i = 0; i < result.alignmentLength; i++)
        {
            if (result.alignment[i] == 0)
            {
                num_match++;
            }
            else if (result.alignment[i] == 3)
            {
                num_mismatch++;
            }
        }
        return (float)num_match / (num_match + num_mismatch);
    }

    inline bool another_chance(const char *target, std::uint16_t t_len, const char *query, std::uint16_t q_len)
    {
        auto result = get_edlib_result(query, q_len, target, t_len, EDLIB_MODE_INFIX, EDLIB_TASK_DISTANCE);
        if ((float)result.editDistance / q_len < 0.15)
        {
            return true;
        }
        return false;
    }

    inline std::string InflateDataAccordingly(std::unique_ptr<biosoup::NucleicAcid> &seq, std::uint32_t start,
                                              std::uint32_t len, bool strand)
    {
        if (!strand)
        {
            return seq->InflateDataReverse(start, len);
        }
        else
        {
            return seq->InflateData(start, len);
        }
    }
    // align the rhs to the lhs
    clipped_alignment<EdlibAlignResult> align_overlap(biosoup::Overlap o, align_reads::Inputs *inputs, const char *target_string, std::uint32_t target_len, std::uint8_t query_group)
    {
        clipped_alignment<EdlibAlignResult> to_return;
        // Get the overlapping segment start and end
        std::uint32_t t_begin = o.lhs_begin;
        std::uint32_t t_end = o.lhs_end;
        std::uint32_t q_begin = o.rhs_begin;
        std::uint32_t q_end = o.rhs_end;

        // Get the query sequence
        auto &query = inputs->get_id_in_group(query_group, o.rhs_id);

        // If the coordinates for q are for the reverse-complement,
        // find the right coordinates for the original strand.
        if (!o.strand)
        {
            q_end = query->inflated_len - o.rhs_begin;
            q_begin = query->inflated_len - o.rhs_end;
        }

        // How much does q extend beyond the ends of t
        // when the overlapping segments are aligned.
        int protrude_left = q_begin - t_begin;
        int protrude_right = (query->inflated_len - q_end) - (target_len - t_end);

        // How much (approximately) to clip off each sequence so that only the overlapping segments
        // are left.
        std::uint32_t q_clip_left;
        std::uint32_t q_clip_right;
        std::uint32_t t_clip_left;
        std::uint32_t t_clip_right;
        if (protrude_left > 0)
        {
            q_clip_left = protrude_left;
            t_clip_left = 0;
        }
        else
        {
            q_clip_left = 0;
            t_clip_left = -protrude_left;
        }
        if (protrude_right > 0)
        {
            q_clip_right = protrude_right;
            t_clip_right = 0;
        }
        else
        {
            q_clip_right = 0;
            t_clip_right = -protrude_right;
        }

        std::uint32_t left_overhang = q_begin - q_clip_left;
        std::uint32_t right_overhang = query->inflated_len - q_clip_right - q_end;
        std::string query_segment;
        auto q_len = query->inflated_len - q_clip_left - q_clip_right;
        auto t_len = target_len - t_clip_left - t_clip_right;
        if (left_overhang > 500)
        {
            query_segment = InflateDataAccordingly(query, q_clip_left + 100, 100, o.strand);
            if (!another_chance(target_string + t_clip_left, 300, query_segment.c_str(), 100))
            {
                to_return.valid = false;
                return to_return;
            }
        }
        if (right_overhang > 500)
        {
            query_segment = InflateDataAccordingly(query, q_clip_left + q_len - 200, 100, o.strand);
            if (!another_chance(target_string + t_clip_left + t_len - 300, 300, query_segment.c_str(), 100))
            {
                to_return.valid = false;
                return to_return;
            }
        }

        query_segment = InflateDataAccordingly(query, q_clip_left, q_len, o.strand);
        to_return.q_start = q_clip_left;
        to_return.q_end = q_clip_left + q_len - 1;
        to_return.t_start = t_clip_left;
        to_return.t_end = t_clip_left + t_len - 1;
        to_return.rhs_id = o.rhs_id;
        to_return.result = get_edlib_result_(query_segment.c_str(), q_len,
                                             target_string + t_clip_left, t_len,
                                             EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
        to_return.identity_score = iden_from_align_result(to_return.result);
        to_return.clipped_query = std::move(query_segment);
        return to_return;
    }

    std::vector<clipped_alignment<EdlibAlignResult>> align_overlaps(std::vector<biosoup::Overlap> &overlaps, std::uint16_t num, align_reads::Inputs &inputs, std::string &target, std::shared_ptr<thread_pool::ThreadPool> &pool, std::uint8_t query_group)
    {
        if (pool == nullptr)
        {
            std::vector<clipped_alignment<EdlibAlignResult>> output;
            output.reserve(num);
            for (auto i = 0; i < num; i++)
            {
                if (overlaps.empty())
                    break;
                auto o = overlaps.back();
                overlaps.pop_back();
                output.push_back(align_overlap(o, &inputs, target.c_str(), target.size(), query_group));
            }
            return output;
        }
        else
        {
            Futures<clipped_alignment<EdlibAlignResult>> futures{pool, num};
            for (auto i = 0; i < num; i++)
            {
                if (overlaps.empty())
                    break;
                auto o = overlaps.back();
                overlaps.pop_back();
                futures.add_inputs(align_overlap, o, &inputs, target.c_str(), target.size(), query_group);
            }
            return futures.get();
        }
    }

    AlignmentSegment get_alignment_segment2(clipped_alignment<EdlibAlignResult> &result, std::string &target_string, bool free_result)
    {
        return {result.clipped_query, 0, target_string, result.t_start, result.result, free_result};
    }

} // namespace align_reads