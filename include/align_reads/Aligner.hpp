#ifndef ALIGN_READS_GENERATOR_HPP_
#define ALIGN_READS_GENERATOR_HPP_

#include "edlib.h"
#include "thread_pool/thread_pool.hpp"

#define EDLIB_MODE_GLOBAL EDLIB_MODE_NW
#define EDLIB_MODE_PREFIX EDLIB_MODE_SHW
#define EDLIB_MODE_INFIX EDLIB_MODE_HW

/*
 * Takes care of alignment operations.
 */

namespace align_reads
{
    // A single alignment task for edlib, from a query to a target sequence
    struct EdlibTask
    {
        const char *query_start;
        const char *target_start;
        std::uint32_t query_len;
        std::uint32_t target_len;
        EdlibAlignMode mode;
        EdlibAlignTask task;

        EdlibTask(const char *q_start, const char *t_start, std::uint32_t q_len, std::uint32_t t_len,
                   EdlibAlignMode mode, EdlibAlignTask task)
            : query_start(q_start), target_start(t_start), query_len(q_len), target_len(t_len), mode(mode), task(task) {}
    };

    /* Stores alignment info of a query w.r.t. to a target.
     * Always exist w.r.t. to a target.
     */
    class AlignmentSegment
    {
    public:
        // q_start/t_start are offsets from the start of query/target
        // when aligning with edlib. Result should contain path.
        // Frees(consumes) the edlib result.
        AlignmentSegment(std::string &query, std::string &target,
                          std::uint32_t q_start, std::uint32_t t_start, EdlibAlignResult &result);
        
        AlignmentSegment(std::string &query, std::string &target,
                          std::uint32_t q_start, std::uint32_t t_start,
                          std::uint32_t q_len, std::uint32_t t_len,
                           EdlibAlignMode mode, EdlibAlignTask task);

        // Get the characters (possibily '_' ) aligned to the target segment.
        std::string &get_aligned_chars();

        // Get the insertion segment at the specified index of the target.
        // -1 for ins segment before 0.
        std::string &get_ins_segment_at(int index);

    private:
        std::uint32_t start_on_target;     // The start index of the aligned segment on the target, 0-based
        std::uint32_t end_on_target; // The end index, 0-based
        std::string aligned_chars;
        std::vector<std::string> ins_segments;

        friend class MultiAlignment;
        friend class AlignmentConverter;
    };

    // Stores info on how a set of sequences align to a single target.
    class MultiAlignment
    {
    public:
        MultiAlignment(std::string &target, std::vector<AlignmentSegment>& segments);
        // steals content of inputs
        MultiAlignment(std::string &&target, std::vector<AlignmentSegment>&& segments);


    private:
        std::string target; // The target sequence
        std::vector<AlignmentSegment> alignment_segments;

        friend class AlignmentConverter;
    };

    // Process multiple edlib tasks, possibly in parallel
    std::vector<EdlibAlignResult> get_edlib_results(std::vector<EdlibTask> &tasks, std::shared_ptr<thread_pool::ThreadPool> &pool);

    // process one edlib task
    EdlibAlignResult get_edlib_result(const char *q_start, const char *t_start,
                                      std::uint32_t q_len, std::uint32_t t_len, EdlibAlignMode mode, EdlibAlignTask task);

    // Integrates alignment, straightaway gives alignment_segment
    AlignmentSegment get_alignment_segment(std::string& query, std::string& target,
         std::uint32_t q_start, std::uint32_t t_start, std::uint32_t q_len, std::uint32_t t_len,
         EdlibAlignMode mode, EdlibAlignTask task);

} // namespace align_reads

#endif // ALIGN_READS_GENERATOR_HPP_