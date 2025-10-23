#include "RearrangeCoverage.h"
#include <iostream>
#include <set>

using namespace std;

RearrangeCoverage::RearrangeCoverage(const AlignmentStore& store) : store(store) {}

void RearrangeCoverage::calculate_coverage(const vector<Event>& events, map<string, int>& lib_coverage_matrix) const
{
    // clear output container
    lib_coverage_matrix.clear();
    
    cout << "calculating coverage for " << events.size() << " events..." << endl;
    
    for (const Event& event : events) {
        lib_coverage_matrix[event.event_id] = calculate_event_coverage(event);
    }
}

int RearrangeCoverage::calculate_event_coverage(const Event& event) const
{
    // skip if contig doesn't exist in this store
    if (!store.has_contig_id(event.contig_id)) {
        return 0;
    }
    
    // create two intervals: [clip_out-1, clip_out] and [clip_in, clip_in+1]
    uint32_t clip_out_0based = event.out_clip - 1; // convert to 0-based
    uint32_t clip_in_0based = event.in_clip - 1;   // convert to 0-based
    
    Interval interval1(event.contig_id, 
                      clip_out_0based > 0 ? clip_out_0based - 1 : 0, 
                      clip_out_0based);
    Interval interval2(event.contig_id, clip_in_0based + 1, clip_in_0based + 2);
    
    // get alignments from both intervals
    auto alignments1 = store.get_alignments_intersecting_interval(interval1);
    auto alignments2 = store.get_alignments_intersecting_interval(interval2);
    
    // collect unique read IDs that pass mutation filter
    set<string> reads1;
    set<string> reads2;
    
    // hardcoded mutation threshold (could be parameterized later)
    const double max_anchor_mutations_percent = 1.0;
    
    // process alignments from first interval
    for (const auto& alignment_ref : alignments1) {
        const Alignment& aln = alignment_ref.get();
        
        // apply mutation filter
        double mutation_percent = (static_cast<double>(aln.get_mutation_count()) / 
            static_cast<double>(aln.read_end - aln.read_start)) * 100.0;
        if (mutation_percent <= max_anchor_mutations_percent) {
            string read_id = store.get_read_id(aln.read_index);
            reads1.insert(read_id);
        }
    }
    
    // process alignments from second interval
    for (const auto& alignment_ref : alignments2) {
        const Alignment& aln = alignment_ref.get();
        
        // apply mutation filter
        double mutation_percent = (static_cast<double>(aln.get_mutation_count()) / 
            static_cast<double>(aln.read_end - aln.read_start)) * 100.0;
        if (mutation_percent <= max_anchor_mutations_percent) {
            string read_id = store.get_read_id(aln.read_index);
            reads2.insert(read_id);
        }
    }
    
    // calculate average of reads covering each interval
    int avg_reads = static_cast<int>((reads1.size() + reads2.size()) / 2);
    return avg_reads;
}
