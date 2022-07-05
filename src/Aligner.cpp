#include <assert.h>

#include "align_reads/Converter.hpp"
#include "edlib.h"
#include "thread_pool/thread_pool.hpp"

#include "align_reads/Aligner.hpp"
#include "align_reads/Utilities.hpp"

namespace align_reads
{
    //------------------inlines----------------------


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

    EdlibAlignResult align_overlap(biosoup::Overlap &o, align_reads::Inputs &inputs, std::string &target)
    {
        // Get the overlapping segment start and end
        std::uint32_t t_begin = o.lhs_begin;
        std::uint32_t t_end = o.lhs_end;
        std::uint32_t q_begin = o.rhs_begin;
        std::uint32_t q_end = o.rhs_end;

        // Get the query sequence
        auto& query = inputs.get_id_in_group(READS_INDEX, o.rhs_id);

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
        int protrude_right = (query->inflated_len - q_end) - (target.size() - t_end);

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
        
        std::string query_segment = query->InflateData(q_clip_left, query->inflated_len - q_clip_left - q_clip_right);
        return get_edlib_result_(query_segment.c_str(), query_segment.size(),
                                        target.c_str() + t_clip_left, target.size() - t_clip_left - t_clip_right,
                                        EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    }

    /*
     std::vector<EdlibAlignResult> align_overlaps(std::vector<biosoup::Overlap>& overlaps, std::uint16_t num)
     {
        Futures<EdlibAlignResult> futures;
        for (auto i = 0; i < num; i++)
        {
            if (overlaps.empty()) break;
            auto& o = overlaps.back();

        }

     } */

} // namespace align_reads