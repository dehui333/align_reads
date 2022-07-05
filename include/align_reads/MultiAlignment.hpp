#ifndef ALIGN_READS_MULTIALIGNMENT_HPP_
#define ALIGN_READS_MULTIALIGNMENT_HPP_

#include "align_reads/AlignmentSegment.hpp"

namespace align_reads
{

    // Stores info on how a set of sequences align to a single target.
    class MultiAlignment
    {
    public:
        MultiAlignment(std::string &target, std::vector<AlignmentSegment> &segments);
        // steals content of inputs
        MultiAlignment(std::string &&target, std::vector<AlignmentSegment> &&segments);

        MultiAlignment(std::string &&target, std::vector<AlignmentSegment> &&segments, std::vector<AlignmentSegment> &&truth);

    private:
        std::string target; // The target sequence
        std::vector<AlignmentSegment> alignment_segments;
        std::vector<AlignmentSegment> truth_to_target;

        friend class AlignmentConverter;
    };

} // namespace align_reads

#endif // ALIGN_READS_MULTIALIGNMENT_HPP_
