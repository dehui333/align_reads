#include <assert.h>

#include "align_reads/Converter.hpp"
#include "edlib.h"
#include "thread_pool/thread_pool.hpp"

#include "align_reads/Aligner.hpp"

#define GAP_CHAR '_'

namespace align_reads
{
    //-----------------AlignmentSegment---------------------------

    AlignmentSegment::AlignmentSegment(std::string &query, std::string &target,
                                         std::uint32_t q_start, std::uint32_t t_start,
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

    //----------------- MultiAlignment---------------------------

    MultiAlignment::MultiAlignment(std::string &target,
     std::vector<AlignmentSegment>& segments) : target(target), alignment_segments(segments) {}
    MultiAlignment::MultiAlignment(std::string &&target,
     std::vector<AlignmentSegment>&& segments) : target(std::move(target)), alignment_segments(std::move(segments)) {}

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

    EdlibAlignResult get_edlib_result(const char *q_start, const char *t_start,
                                      std::uint32_t q_len, std::uint32_t t_len, EdlibAlignMode mode, EdlibAlignTask task)
    {
        return edlibAlign(q_start, q_len, t_start, t_len, edlibNewAlignConfig(-1, mode, task, NULL, 0));
    }

    AlignmentSegment get_alignment_segment(std::string &query, std::string &target,
                                            std::uint32_t q_start, std::uint32_t t_start, std::uint32_t q_len, std::uint32_t t_len,
                                            EdlibAlignMode mode, EdlibAlignTask task)
    {
        auto result = edlibAlign(query.c_str() + q_start, q_len,
                                 target.c_str() + t_start, t_len, edlibNewAlignConfig(-1, mode, task, NULL, 0));
        return {query, target, q_start, t_start, result};
    }

} // namespace align_reads