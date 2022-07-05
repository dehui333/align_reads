
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
        : target(std::move(target)), alignment_segments(std::move(segments)), truth_to_target(std::move(truth)) {}

} // namespace align_reads