#ifndef ALIGN_READS_GENERATOR_HPP_
#define ALIGN_READS_GENERATOR_HPP_

#include "biosoup/overlap.hpp"
#include "edlib.h"
#include "thread_pool/thread_pool.hpp"

#include "align_reads/AlignmentSegment.hpp"
#include "align_reads/Inputs.hpp"
#include "align_reads/MultiAlignment.hpp"

#define EDLIB_MODE_GLOBAL EDLIB_MODE_NW
#define EDLIB_MODE_PREFIX EDLIB_MODE_SHW
#define EDLIB_MODE_INFIX EDLIB_MODE_HW




// !!!!!!!!! biosoup::Nucleic_acid's ReverseAndComplement() changes the
// state of the sequence object and is not thread safe

// todo?? : local realignment to have affine gap...?
// todo?? : split alignment into subwindows by using estimation from overlap...?


/*
 * Takes care of alignment operations.
 *
 *  - Gives friend access to Converter.
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

        EdlibTask(const char *q_start, std::uint32_t q_len, const char *t_start, std::uint32_t t_len,
                  EdlibAlignMode mode, EdlibAlignTask task)
            : query_start(q_start), target_start(t_start), query_len(q_len), target_len(t_len), mode(mode), task(task) {}
    };

    // Gives info on which segments of the query and target are aligned
    template <typename R>
    struct clipped_alignment
    {
        std::string clipped_query;
        std::uint32_t q_start;
        std::uint32_t q_end; // ends are inclusive
        std::uint32_t t_start;
        std::uint32_t t_end;
        R result;

    };


    // Process multiple edlib tasks, possibly in parallel
    std::vector<EdlibAlignResult> get_edlib_results(std::vector<EdlibTask> &tasks, std::shared_ptr<thread_pool::ThreadPool> &pool);

    // process one edlib task
    EdlibAlignResult get_edlib_result(const char *q_start, std::uint32_t q_len,
                                      const char *t_start, std::uint32_t t_len, EdlibAlignMode mode, EdlibAlignTask task);

    // align some number of overlapping segments.
    // will remove the aligned ones from the container.
    std::vector<clipped_alignment<EdlibAlignResult>> align_overlaps(std::vector<biosoup::Overlap> &overlaps, std::uint16_t num, align_reads::Inputs &inputs, std::string &target, std::uint8_t query_group);

    // Integrates alignment, straightaway gives alignment_segment
    AlignmentSegment get_alignment_segment(std::string &query, std::uint32_t q_start, std::uint32_t q_len,
                                           std::string &target, std::uint32_t t_start, std::uint32_t t_len,
                                           EdlibAlignMode mode, EdlibAlignTask task);

} // namespace align_reads

#endif // ALIGN_READS_GENERATOR_HPP_