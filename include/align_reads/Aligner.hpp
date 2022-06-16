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
    struct edlib_task
    {
        const char *query_start;
        const char *target_start;
        std::uint32_t query_len;
        std::uint32_t target_len;
        EdlibAlignMode mode;
        EdlibAlignTask task;

        edlib_task(const char *q_start, const char *t_start, std::uint32_t q_len, std::uint32_t t_len,
                   EdlibAlignMode mode, EdlibAlignTask task)
            : query_start(q_start), target_start(t_start), query_len(q_len), target_len(t_len), mode(mode), task(task) {}
    };

    // Process multiple edlib tasks, possibly in parallel
    std::vector<EdlibAlignResult> get_edlib_results(std::vector<edlib_task> &tasks, std::shared_ptr<thread_pool::ThreadPool> &pool);

    // process one edlib task
    EdlibAlignResult get_edlib_result(const char *q_start, const char *t_start,
                                      std::uint32_t q_len, std::uint32_t t_len, EdlibAlignMode mode, EdlibAlignTask task);

} // namespace align_reads

#endif // ALIGN_READS_GENERATOR_HPP_