#ifndef REARRANGE_COVERAGE_H
#define REARRANGE_COVERAGE_H

#include "alignment_store.h"
#include "RearrangeEventGroup.h"
#include <map>
#include <string>

// worker class for calculating coverage for events
class RearrangeCoverage {
private:
    const AlignmentStore& store;

public:
    RearrangeCoverage(const AlignmentStore& store);
    
    // calculate coverage for a set of events - fills output container by reference
    void calculate_coverage(const std::vector<Event>& events, std::map<std::string, int>& lib_coverage_matrix) const;

private:
    // calculate coverage for a single event
    int calculate_event_coverage(const Event& event) const;
};

#endif // REARRANGE_COVERAGE_H
