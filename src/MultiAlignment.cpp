#include "align_reads/MultiAlignment.hpp"

namespace align_reads
{

    //----------------- MultiAlignment---------------------------

    MultiAlignment::MultiAlignment(std::string &target,
                                   std::vector<AlignmentSegment> &segments) : target(target), alignment_segments(segments) {}
    MultiAlignment::MultiAlignment(std::string &&target,
                                   std::vector<AlignmentSegment> &&segments) : target(std::move(target)), alignment_segments(std::move(segments)) {}
    MultiAlignment::MultiAlignment(std::string &&target,
                                   std::vector<AlignmentSegment> &&segments,
                                   std::vector<AlignmentSegment> &&truth)
        : target(std::move(target)), alignment_segments(std::move(segments)), truth_to_target(std::move(truth)) 
        {
            largest_truth_start = 0;
            smallest_truth_end = -1;
            for (auto& s : truth_to_target)
            {
                if (s.start_on_target > largest_truth_start) largest_truth_start = s.start_on_target;
                if (s.end_on_target < smallest_truth_end) smallest_truth_end = s.end_on_target;
            }
        }

} // namespace align_reads