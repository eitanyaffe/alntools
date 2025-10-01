#include "RearrangeReadEvent.h"
#include "RearrangeVerify.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <map>

using namespace std;

// convert EventTestResult to string for rejection counting
string event_test_result_to_string(EventTestResult result)
{
    switch (result) {
    case EventTestResult::FOUND_EVENT:
        return "found_event";
    case EventTestResult::REJECT_ALIGNMENT_OVERLAP:
        return "reject_alignment_overlap";
    case EventTestResult::REJECT_INSERT_ELEMENT_REINSERTED:
        return "reject_insert_element_reinserted";
    case EventTestResult::REJECT_INSERT_ELEMENT_OVERLAPS:
        return "reject_insert_element_overlaps";
    case EventTestResult::REJECT_ELEMENT_TOO_SHORT:
        return "reject_element_too_short";
    case EventTestResult::REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE:
        return "reject_contained_element_margin_too_large";
    case EventTestResult::REJECT_READ_GAP_TOO_LARGE:
        return "reject_read_gap_too_large";
    case EventTestResult::REJECT_INSERT_GAP_TOO_LARGE:
        return "reject_insert_gap_too_large";
    default:
        return "unknown";
    }
}

RearrangeReadEvent::RearrangeReadEvent(const AlignmentStore& store,
                                     RearrangeVerify* verifier,
                                     ResolveSeams resolve_seams,
                                     const AssemblySequences* assembly_seqs,
                                     const ReadSequences* read_seqs,
                                     int max_margin,
                                     int min_element_length,
                                     int min_anchor_length,
                                     double max_anchor_mutations_percent,
                                     double max_element_mutation_percent)
    : store(store)
    , max_margin(max_margin)
    , min_element_length(min_element_length)
    , min_anchor_length(min_anchor_length)
    , max_anchor_mutations_percent(max_anchor_mutations_percent)
    , max_element_mutation_percent(max_element_mutation_percent)
    , verifier(verifier)
    , resolve_seams(resolve_seams)
    , assembly_sequences(assembly_seqs)
    , read_sequences(read_seqs)
{
}

void RearrangeReadEvent::detect_events(const string& lib_id, vector<ReadEvent>& output_events, 
                                      map<string, size_t>& output_rejection_counts, const AlignmentStore& store,
                                      const string& output_prefix)
{
    cout << "detecting rearrangement events for library " << lib_id << "..." << endl;
    current_lib_id = lib_id;
    
    // debugging settings
    const string focal_read_id = "XXX_m84221_250708_233611_s3/168956412/ccs";
    
    const auto& reads = store.get_reads();
    const auto& alignments = store.get_alignments();
    
    // clear output containers
    output_events.clear();
    output_rejection_counts.clear();
    
    size_t tested_reads = 0;
    size_t reads_with_events = 0;
    
    // determine which reads to process
    set<uint32_t> reads_to_process;
    if (relevant_read_indices.empty()) {
        // no interval filtering - process all reads
        for (size_t read_idx = 0; read_idx < reads.size(); ++read_idx) {
            reads_to_process.insert(static_cast<uint32_t>(read_idx));
        }
        cout << "processing " << reads.size() << " reads (no interval filtering)..." << endl;
    } else {
        // interval filtering - process only relevant reads
        reads_to_process = relevant_read_indices;
        cout << "processing " << reads_to_process.size() << " reads (interval filtered)..." << endl;
    }
    
    
    // process relevant reads
    for (uint32_t read_idx : reads_to_process) {
        
        // get alignments for this read (already sorted by read_start)
        auto it = store.get_read_to_alignment_indices().find(static_cast<uint32_t>(read_idx));
        if (it == store.get_read_to_alignment_indices().end()) {
            continue; // no alignments for this read
        }
        
        const vector<size_t>& alignment_indices = it->second;
        
        // skip reads with less than 2 alignments
        if (alignment_indices.size() < 2) {
            continue;
        }
        
        // check if this is the focal read for verbose debugging
        string current_read_id = store.get_read_id(read_idx);
        bool verbose = !focal_read_id.empty() && current_read_id == focal_read_id;
        
        if (verbose) {
            cout << "VERBOSE: starting to process focal read: " << current_read_id << endl;
        }
        
        // process this read
        if (process_read(alignment_indices, alignments, output_events, output_rejection_counts, verbose)) {
            reads_with_events++;
        }
        
        tested_reads++;
    }
    
    cout << "tested " << tested_reads << " reads" << endl;
    cout << "found " << output_events.size() << " events in " << reads_with_events << " reads" << endl;
    
    // get read sequences if needed for seam resolution or verification
    unordered_map<string, string> sequences;
    bool need_sequences = (read_sequences && (resolve_seams == ResolveSeams::READS_ONLY || 
                                            resolve_seams == ResolveSeams::COMPLETE)) || verifier;
    
    if (need_sequences && !output_events.empty()) {
        get_read_sequences(lib_id, output_events, sequences, output_prefix);
    }
    
    // resolve read seams at the end
    if (read_sequences && (resolve_seams == ResolveSeams::READS_ONLY || 
                          resolve_seams == ResolveSeams::COMPLETE)) {
        resolve_read_seams(lib_id, output_events, sequences);
    }
    
    // verify all detected events at the end
    if (verifier) {
        verifier->verify_all_events(lib_id, output_events, store, sequences);
    }
    
    // print rejection summary
    size_t total_rejected = 0;
    for (const auto& pair : output_rejection_counts) {
        total_rejected += pair.second;
        cout << "rejected " << pair.second << " candidates due to: " << pair.first << endl;
    }
    cout << "total rejected: " << total_rejected << " candidates" << endl;
}

void RearrangeReadEvent::set_intervals(const vector<Interval>& intervals)
{
    this->intervals = intervals;
    relevant_read_indices.clear();
    relevant_alignments.clear();
    
    if (intervals.empty()) {
        cout << "no intervals provided - will process all reads" << endl;
        return;
    }
    
    cout << "applying interval filtering for " << intervals.size() << " intervals..." << endl;
    
    // find all alignments that overlap with any interval
    for (const Interval& interval : intervals) {
        auto alignments_in_interval = store.get_alignments_intersecting_interval(interval);
        
        for (const auto& alignment_ref : alignments_in_interval) {
            const Alignment& aln = alignment_ref.get();
            relevant_alignments.insert(&aln);
            relevant_read_indices.insert(aln.read_index);
        }
    }
    
    cout << "found " << relevant_alignments.size() << " relevant alignments from " 
         << relevant_read_indices.size() << " reads" << endl;
}

size_t RearrangeReadEvent::get_total_rejected_count(const map<string, size_t>& rejection_counts) const
{
    size_t total = 0;
    for (const auto& pair : rejection_counts) {
        total += pair.second;
    }
    return total;
}

void RearrangeReadEvent::print_event(const ReadEvent& event, const AlignmentStore& store)
{
    cout << "EVENT: read_id=" << event.read_id 
         << " type=" << event.type
         << " contig=" << event.contig_id 
         << " out_clip=" << event.out_clip 
         << " in_clip=" << event.in_clip
         << " read_strand=" << event.read_strand;
    
    if (!event.element_contig.empty()) {
        cout << " element_contig=" << event.element_contig
             << " element_strand=" << event.element_strand
             << " element_start=" << event.element_start
             << " element_end=" << event.element_end;
        
        if (!event.read_seams.empty()) {
            cout << " read_seams_len=" << event.read_seams.length();
        }
        if (!event.assembly_seams.empty()) {
            cout << " assembly_seams_len=" << event.assembly_seams.length();
        }
    }
    
    cout << endl;
    
    // print alignment details using print_alignment
    print_alignment(event.alignment_A, store, "  Alignment A");
    print_alignment(event.alignment_B, store, "  Alignment B");
    
    if (event.has_alignment_X) {
        print_alignment(event.alignment_X, store, "  Alignment X");
    }
}

void RearrangeReadEvent::print_read_alignments(const vector<size_t>& alignment_indices,
                                              const vector<Alignment>& alignments,
                                              const vector<bool>& is_valid_anchor_vec) const
{
    cout << "VERBOSE: read has " << alignment_indices.size() << " alignments:" << endl;
    for (size_t i = 0; i < alignment_indices.size(); ++i) {
        const Alignment& aln = alignments[alignment_indices[i]];
        cout << "  [" << i << "] read:" << aln.read_start << "-" << aln.read_end 
             << " contig:" << store.get_contig_id(aln.contig_index) << ":" 
             << aln.contig_start << "-" << aln.contig_end
             << " strand:" << (aln.is_reverse ? "-" : "+")
             << " valid_anchor:" << (is_valid_anchor_vec[i] ? "yes" : "no") << endl;
    }
}

void RearrangeReadEvent::add_rejection(EventTestResult result, map<string, size_t>& rejection_counts, bool verbose)
{
    string reason = event_test_result_to_string(result);
    rejection_counts[reason]++;
    if (verbose) {
        cout << "rejecting: " << reason << endl;
    }
}

bool RearrangeReadEvent::is_valid_anchor(const Alignment& aln, bool verbose, size_t index) const
{
    // check minimum anchor length
    uint32_t alignment_length = aln.read_end - aln.read_start;
    if (static_cast<int>(alignment_length) < min_anchor_length) {
        if (verbose) {
            cout << "VERBOSE: alignment [" << index << "] rejected as anchor - too short (length=" << alignment_length << ", min=" << min_anchor_length << ")" << endl;
        }
        return false;
    }
    
    // check mutation rate
    double mutation_percent = static_cast<double>(aln.get_mutation_count()) / 
        static_cast<double>(alignment_length);
    if (mutation_percent > max_anchor_mutations_percent) {
        if (verbose) {
            cout << "VERBOSE: alignment [" << index << "] rejected as anchor - too many mutations (rate=" << mutation_percent << ", max=" << max_anchor_mutations_percent << ")" << endl;
        }
        return false;
    }
    
    // check alignment overlap
    uint32_t overlap = store.get_alignment_overlap(aln);
    if (overlap > static_cast<uint32_t>(max_margin)) {
        if (verbose) {
            cout << "VERBOSE: alignment [" << index << "] rejected as anchor - overlap too large (overlap=" << overlap << ", max=" << max_margin << ")" << endl;
        }
        return false;
    }
    
    if (verbose) {
        cout << "VERBOSE: alignment [" << index << "] accepted as valid anchor" << endl;
    }
    
    return true;
}

bool RearrangeReadEvent::is_valid_element(const Alignment& aln) const
{
    // check minimum element length
    uint32_t alignment_length = aln.read_end - aln.read_start;
    if (static_cast<int>(alignment_length) < min_element_length) {
        return false;
    }
    
    // check mutation rate
    double mutation_percent = static_cast<double>(aln.get_mutation_count()) / 
        static_cast<double>(alignment_length);
    if (mutation_percent > max_element_mutation_percent) {
        return false;
    }
    
    return true;
}

bool RearrangeReadEvent::process_read(const vector<size_t>& alignment_indices, 
                                     const vector<Alignment>& alignments,
                                     vector<ReadEvent>& output_events,
                                     map<string, size_t>& rejection_counts,
                                     bool verbose)
{
    bool found_event_in_read = false;
    
    // pre-processing: calculate valid anchors for all alignments
    vector<bool> is_valid_anchor_vec(alignment_indices.size(), false);
    for (size_t i = 0; i < alignment_indices.size(); ++i) {
        const Alignment& aln = alignments[alignment_indices[i]];
        is_valid_anchor_vec[i] = is_valid_anchor(aln, verbose, i);
    }
    
    // print all alignments if verbose
    if (verbose) {
        print_read_alignments(alignment_indices, alignments, is_valid_anchor_vec);
    }
    
    
    // try trio events first (A-X-B pattern)
    if (process_read_trios(alignment_indices, alignments, is_valid_anchor_vec, 
                          output_events, rejection_counts, verbose)) {
        found_event_in_read = true;
    }
    
    // try pair events (A-B pattern with gap)
    if (process_read_pairs(alignment_indices, alignments, is_valid_anchor_vec,
                          output_events, rejection_counts, verbose)) {
        found_event_in_read = true;
    }
    
    return found_event_in_read;
}

bool RearrangeReadEvent::process_read_trios(const vector<size_t>& alignment_indices,
                                           const vector<Alignment>& alignments,
                                           const vector<bool>& is_valid_anchor_vec,
                                           vector<ReadEvent>& output_events,
                                           map<string, size_t>& rejection_counts,
                                           bool verbose)
{
    bool found_event = false;

    // test trio events - based on original sequential logic
    if (verbose) {
        cout << "VERBOSE: starting trio event loop" << endl;
    }
    
    for (size_t i = 0; i < alignment_indices.size(); ++i) {
        // anchor1: current alignment, must be valid anchor
        if (!is_valid_anchor_vec[i]) {
            if (verbose) {
                cout << "VERBOSE: trio i=" << i << " - skipping, not valid anchor" << endl;
            }
            continue;
        }
        const Alignment& anchor1 = alignments[alignment_indices[i]];
        
        // early reject if anchor1 is not in relevant alignments (when interval filtering is active)
        if (!relevant_alignments.empty() && relevant_alignments.find(&anchor1) == relevant_alignments.end()) {
            if (verbose) {
                cout << "VERBOSE: trio i=" << i << " - skipping, anchor1 not in relevant alignments" << endl;
            }
            continue;
        }
        
        // find anchor2: first valid anchor after anchor1, skipping the very next alignment
        size_t j = i + 2; // skip i+1
        while (j < alignment_indices.size() && !is_valid_anchor_vec[j]) {
            j++;
        }
        if (j >= alignment_indices.size()) {
            if (verbose) {
                cout << "VERBOSE: trio i=" << i << " - skipping, no valid anchor2 found" << endl;
            }
            continue;
        }
        const Alignment& anchor2 = alignments[alignment_indices[j]];
        
        // early reject if anchor2 is not in relevant alignments (when interval filtering is active)
        if (!relevant_alignments.empty() && relevant_alignments.find(&anchor2) == relevant_alignments.end()) {
            if (verbose) {
                cout << "VERBOSE: trio i=" << i << ", j=" << j << " - skipping, anchor2 not in relevant alignments" << endl;
            }
            continue;
        }
                 
        // verify that anchors are on same strand and contig
        if (anchor1.is_reverse != anchor2.is_reverse || anchor1.contig_index != anchor2.contig_index) {
            if (verbose) {
                cout << "VERBOSE: trio i=" << i << ", j=" << j << " - skipping, different strand/contig" << endl;
            }
            continue;
        }
        
        // find element: best alignment in interval between anchor1 and anchor2
        uint32_t interval_start = anchor1.read_end;
        uint32_t interval_end = anchor2.read_start;
        
        const Alignment* element = store.get_alignment_in_read(anchor1.read_index, interval_start, interval_end);
        if (!element || !is_valid_element(*element)) {
            if (verbose) {
                cout << "VERBOSE: trio i=" << i << ", j=" << j << " - skipping, no valid element found" << endl;
            }
            continue;
        }
        
        // early reject if element is not in relevant alignments (when interval filtering is active)
        if (!relevant_alignments.empty() && relevant_alignments.find(element) == relevant_alignments.end()) {
            if (verbose) {
                cout << "VERBOSE: trio i=" << i << ", j=" << j << " - skipping, element not in relevant alignments" << endl;
            }
            continue;
        }
        
        // test gaps (allowing max_margin overlap/gap with anchor1 and anchor2)
        int32_t delta1 = static_cast<int32_t>(element->read_start) - static_cast<int32_t>(anchor1.read_end);
        int32_t delta2 = static_cast<int32_t>(anchor2.read_start) - static_cast<int32_t>(element->read_end);
        
        if (abs(delta1) > max_margin || abs(delta2) > max_margin) {
            if (verbose) {
                cout << "VERBOSE: trio i=" << i << ", j=" << j << " - skipping, gap too large (delta1=" << delta1 << ", delta2=" << delta2 << ")" << endl;
            }
            continue;
        }
        
        // define A, B, X (after all other validation)
        bool is_read_reversed = anchor1.is_reverse;  // determine flag from anchor1/anchor2 strand
        const Alignment* A = is_read_reversed ? &anchor2 : &anchor1;
        const Alignment* B = is_read_reversed ? &anchor1 : &anchor2;
        
        // create a copy of element to potentially modify its strand
        Alignment X = *element;
        
        // adjust X strand relative to read orientation if needed
        if (is_read_reversed) {
            X.is_reverse = !X.is_reverse;
        }
        
        // test A B order
        if (A->contig_start >= B->contig_start) {
            if (verbose) {
                cout << "VERBOSE: trio i=" << i << ", j=" << j << " - skipping, A not left of B on contig" << endl;
            }
            continue;
        }
        
        if (verbose) {
            cout << "VERBOSE: testing trio event at indices i=" << i << ", j=" << j << endl;
        }
        
        ReadEvent event;
        EventTestResult result = test_trio_event(*A, *B, X, is_read_reversed, event, rejection_counts, verbose);
        
        if (verbose) {
            cout << "VERBOSE: trio event result: " << event_test_result_to_string(result) << endl;
        }
        
        if (result == EventTestResult::FOUND_EVENT) {
            // temporarily skip inline verification - will be moved to end of detect_events
            // if (!verifier || verifier->verify_event(event, current_lib_id, store) == VerifyResult::SUCCESS) {
                output_events.push_back(event);
                found_event = true;
                if (verbose) {
                    cout << "VERBOSE: trio event accepted and added to output" << endl;
                }
            // }
        }
    }
    
    return found_event;
}

bool RearrangeReadEvent::process_read_pairs(const vector<size_t>& alignment_indices,
                                           const vector<Alignment>& alignments,
                                           const vector<bool>& is_valid_anchor_vec,
                                           vector<ReadEvent>& output_events,
                                           map<string, size_t>& rejection_counts,
                                           bool verbose)
{
    bool found_event = false;
    
    // test pair events - consecutive anchors only
    if (verbose) {
        cout << "VERBOSE: starting pair event loop" << endl;
    }
    
    for (size_t i = 0; i + 1 < alignment_indices.size(); ++i) {
        // both alignments must be valid anchors
        if (!is_valid_anchor_vec[i] || !is_valid_anchor_vec[i + 1]) {
            if (verbose) {
                cout << "VERBOSE: pair i=" << i << " - skipping, not both valid anchors" << endl;
            }
            continue;
        }
        
        const Alignment& anchor1 = alignments[alignment_indices[i]];
        const Alignment& anchor2 = alignments[alignment_indices[i + 1]];
        
        // early reject if either anchor is not in relevant alignments (when interval filtering is active)
        if (!relevant_alignments.empty()) {
            if (relevant_alignments.find(&anchor1) == relevant_alignments.end()) {
                if (verbose) {
                    cout << "VERBOSE: pair i=" << i << " - skipping, anchor1 not in relevant alignments" << endl;
                }
                continue;
            }
            if (relevant_alignments.find(&anchor2) == relevant_alignments.end()) {
                if (verbose) {
                    cout << "VERBOSE: pair i=" << i << " - skipping, anchor2 not in relevant alignments" << endl;
                }
                continue;
            }
        }
        
        // test strand and contig
        if (anchor1.is_reverse != anchor2.is_reverse || anchor1.contig_index != anchor2.contig_index) {
            if (verbose) {
                cout << "VERBOSE: pair i=" << i << " - skipping, different strand/contig" << endl;
            }
            continue;
        }
        
        // test gaps (allowing max_margin overlap/gap between anchors)
        int32_t delta = static_cast<int32_t>(anchor2.read_start) - static_cast<int32_t>(anchor1.read_end);
        
        if (abs(delta) > max_margin) {
            if (verbose) {
                cout << "VERBOSE: pair i=" << i << " - skipping, gap too large (delta=" << delta << ")" << endl;
            }
            continue;
        }
        
        // define A, B (after all other validation)
        bool is_read_reversed = anchor1.is_reverse;  // determine flag from anchor1/anchor2 strand
        const Alignment* A = is_read_reversed ? &anchor2 : &anchor1;
        const Alignment* B = is_read_reversed ? &anchor1 : &anchor2;
        
        // test A B order
        if (A->contig_start >= B->contig_start) {
            if (verbose) {
                cout << "VERBOSE: pair i=" << i << " - skipping, A not left of B on contig" << endl;
            }
            continue;
        }
        
        if (verbose) {
            cout << "VERBOSE: testing pair event at index i=" << i << endl;
        }
        
        ReadEvent event;
        EventTestResult result = test_pair_event(*A, *B, is_read_reversed, event, rejection_counts, verbose);
        
        if (verbose) {
            cout << "VERBOSE: pair event result: " << event_test_result_to_string(result) << endl;
        }
        
        if (result == EventTestResult::FOUND_EVENT) {
            // temporarily skip inline verification - will be moved to end of detect_events
            // if (!verifier || verifier->verify_event(event, current_lib_id, store) == VerifyResult::SUCCESS) {
                output_events.push_back(event);
                found_event = true;
                if (verbose) {
                    cout << "VERBOSE: pair event accepted and added to output" << endl;
                }
            // }
        }
    }
    
    return found_event;
}

EventTestResult RearrangeReadEvent::test_trio_event(const Alignment& A, const Alignment& B, const Alignment& X,
                                                   bool is_read_reversed, ReadEvent& event,
                                                   map<string, size_t>& rejection_counts, bool verbose)
{
    // assert A is left of B on contig
    if (A.contig_start >= B.contig_start) {
        cerr << "error: assertion failed - A must be left of B" << endl;
        exit(1);
    }
    
    // assert alignments are properly ordered
    if (A.contig_start >= A.contig_end || B.contig_start >= B.contig_end || X.contig_start >= X.contig_end) {
        cerr << "error: assertion failed - alignment coordinates malformed" << endl;
        exit(1);
    }
    
    // check for alignment overlaps
    if (alignments_overlap(A, B, verbose)) {
        add_rejection(EventTestResult::REJECT_ALIGNMENT_OVERLAP, rejection_counts, verbose);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    if (alignments_overlap(A, X, verbose)) {
        add_rejection(EventTestResult::REJECT_ALIGNMENT_OVERLAP, rejection_counts, verbose);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    if (alignments_overlap(X, B, verbose)) {
        add_rejection(EventTestResult::REJECT_ALIGNMENT_OVERLAP, rejection_counts, verbose);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    // element is outside span → insertion
    if (X.contig_end <= A.contig_end || X.contig_start >= B.contig_start || X.contig_index != A.contig_index) {
        // check if A and B are close enough (A is on the left) for insertion
        uint32_t gap_A_B = B.contig_start - A.contig_end;
        if (gap_A_B > static_cast<uint32_t>(max_margin)) {
            add_rejection(EventTestResult::REJECT_INSERT_GAP_TOO_LARGE, rejection_counts, verbose);
            return EventTestResult::REJECT_INSERT_GAP_TOO_LARGE;
        }
        
        event = create_event(A, B, &X, RearrangementType::LARGE_INSERT, is_read_reversed, verbose);
        return EventTestResult::FOUND_EVENT;
    }
    
    // element is inside span - calculate max distance between X and A, and X and B
    uint32_t distance_X_A = X.contig_start - A.contig_end;
    uint32_t distance_X_B = B.contig_start - X.contig_end;
    uint32_t max_distance = max(distance_X_A, distance_X_B);
    
    if (max_distance > static_cast<uint32_t>(max_margin)) {
        add_rejection(EventTestResult::REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE, rejection_counts, verbose);
        return EventTestResult::REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE;
    }
    
    // check strand orientation
    if (X.is_reverse) {
        // X is on different strand from read orientation → insertion
        event = create_event(A, B, &X, RearrangementType::LARGE_INVERT, is_read_reversed, verbose);
        return EventTestResult::FOUND_EVENT;
    } else {
        // X is on same strand as read orientation → reject as re-insertion
        add_rejection(EventTestResult::REJECT_INSERT_ELEMENT_REINSERTED, rejection_counts, verbose);
        return EventTestResult::REJECT_INSERT_ELEMENT_REINSERTED;
    }
}

EventTestResult RearrangeReadEvent::test_pair_event(const Alignment& A, const Alignment& B,
                                                   bool is_read_reversed, ReadEvent& event,
                                                   map<string, size_t>& rejection_counts, bool verbose)
{
    // assert A is left of B on contig
    if (A.contig_start >= B.contig_start) {
        cerr << "error: assertion failed - A must be left of B" << endl;
        exit(1);
    }
    
    // assert alignments are properly ordered
    if (A.contig_start >= A.contig_end || B.contig_start >= B.contig_end) {
        cerr << "error: assertion failed - alignment coordinates malformed" << endl;
        exit(1);
    }
    
    // check for alignment overlap
    if (alignments_overlap(A, B, verbose)) {
        add_rejection(EventTestResult::REJECT_ALIGNMENT_OVERLAP, rejection_counts, verbose);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    // check gap in read coords between A and B
    uint32_t gap_A_B = (B.read_start > A.read_end) ? (B.read_start - A.read_end) : 0;
    if (gap_A_B > static_cast<uint32_t>(max_margin)) {
        add_rejection(EventTestResult::REJECT_READ_GAP_TOO_LARGE, rejection_counts, verbose);
        return EventTestResult::REJECT_READ_GAP_TOO_LARGE;
    }
    
    // check if inferred element X is long enough
    uint32_t element_length = B.contig_start - A.contig_end;
    if (static_cast<int>(element_length) < min_element_length) {
        add_rejection(EventTestResult::REJECT_ELEMENT_TOO_SHORT, rejection_counts, verbose);
        return EventTestResult::REJECT_ELEMENT_TOO_SHORT;
    }
    
    // deletion event
    event = create_event(A, B, nullptr, RearrangementType::LARGE_DELETE, is_read_reversed, verbose);
    return EventTestResult::FOUND_EVENT;
}


string RearrangeReadEvent::alignment_to_strand_string(bool is_reverse) const
{
    return is_reverse ? "-" : "+";
}

bool RearrangeReadEvent::alignments_overlap(const Alignment& A, const Alignment& B, bool verbose) const
{
    // calculate overlap in read coordinates
    int32_t read_overlap = 0;
    if (A.read_start < B.read_end && B.read_start < A.read_end) {
        uint32_t overlap_start = std::max(A.read_start, B.read_start);
        uint32_t overlap_end = std::min(A.read_end, B.read_end);
        read_overlap = static_cast<int32_t>(overlap_end - overlap_start);
    }
    
    // calculate overlap in contig coordinates (if on same contig)
    int32_t contig_overlap = 0;
    if (A.contig_index == B.contig_index) {
        if (A.contig_start < B.contig_end && B.contig_start < A.contig_end) {
            uint32_t overlap_start = std::max(A.contig_start, B.contig_start);
            uint32_t overlap_end = std::min(A.contig_end, B.contig_end);
            contig_overlap = static_cast<int32_t>(overlap_end - overlap_start);
        }
    }
    
    bool overlaps = (read_overlap > max_margin) || (contig_overlap > max_margin);
    
    if (verbose) {
        cout << "checking alignment overlap:" << endl;
        print_alignment(A, store, "  A");
        print_alignment(B, store, "  B");
        cout << "  read_overlap: " << read_overlap << ", contig_overlap: " << contig_overlap 
             << ", max_margin: " << max_margin << ", overlaps: " << (overlaps ? "true" : "false") << endl;
    }
    
    return overlaps;
}


SeamInterval RearrangeReadEvent::create_seam_interval(const Alignment& align1, const Alignment& align2,
                                                    const string& lib_id, bool is_read_reversed) const
{
    // read coordinates - switch alignments if read is reversed
    uint32_t start, end;
    if (is_read_reversed) {
        massert(align2.read_start < align1.read_end, "align2.read_start must be less than align1.read_end");
        start = align2.read_end;
        end = align1.read_start;
    } else {
        massert(align1.read_start < align2.read_end, "align1.read_start must be less than align2.read_end");
        start = align1.read_end;
        end = align2.read_start;
    }
    
    // validate margin
    int32_t margin = abs(static_cast<int32_t>(end) - static_cast<int32_t>(start));
    if (margin > max_margin) {
        cerr << "read seam exceeds max_margin: " << margin << " > " << max_margin << endl;
        print_alignment(align1, store, "align1");
        print_alignment(align2, store, "align2");
        cerr << "start=" << start << " end=" << end << " is_read_reversed=" << is_read_reversed << endl;
    }
    massert(margin <= max_margin, "read seam exceeds max_margin");
    
    // determine if this is a gap or overlap using same logic as extract_seam
    bool is_gap = (start < end);  // gap if start < end, overlap if start >= end
    
    string read_id = store.get_read_id(align1.read_index);
    return SeamInterval(lib_id, read_id, start, end, is_read_reversed, is_gap);
}

void RearrangeReadEvent::get_read_sequences(const string& lib_id, const vector<ReadEvent>& events, 
                                           unordered_map<string, string>& sequences, const string& output_prefix)
{
    // collect all unique read IDs from events (for both seam resolution and verification)
    unordered_set<string> unique_read_ids;
    
    for (const ReadEvent& event : events) {
        // always collect read ID for verification
        unique_read_ids.insert(event.read_id);
        
        // also collect read IDs from seam intervals if present
        if (!event.read_seam_intervals.empty()) {
            for (const SeamInterval& interval : event.read_seam_intervals) {
                massert(interval.lib_id == lib_id, 
                       "interval lib_id (%s) does not match expected lib_id (%s)", 
                       interval.lib_id.c_str(), lib_id.c_str());
                unique_read_ids.insert(interval.read_id);
            }
        }
    }
    
    if (unique_read_ids.empty()) {
        cout << "no read sequences needed for library " << lib_id << endl;
        return;
    }
    
    cout << "getting sequences for " << unique_read_ids.size() << " unique reads..." << endl;
    
    // convert set to vector for ReadSequences interface
    vector<string> read_ids(unique_read_ids.begin(), unique_read_ids.end());
    
    // get sequences with caching support
    read_sequences->get_sequences_with_cache(lib_id, read_ids, sequences, output_prefix);
}

void RearrangeReadEvent::resolve_read_seams(const string& lib_id, vector<ReadEvent>& events, 
                                           const unordered_map<string, string>& sequences)
{
    cout << "resolving read seams for library " << lib_id << "..." << endl;
    
    size_t events_with_seams = 0;
    size_t resolved_events = 0;
    
    // process each event directly
    for (ReadEvent& event : events) {
        if (!event.read_seam_intervals.empty()) {
            events_with_seams++;
            
            // build read_seams string from intervals
            vector<string> seam_parts;
            for (const SeamInterval& interval : event.read_seam_intervals) {
                auto seq_it = sequences.find(interval.read_id);
                if (seq_it != sequences.end()) {
                    string seam = extract_sequence_interval(seq_it->second, interval);
                    // add prefix based on gap/overlap type
                    if (!seam.empty()) {
                        if (interval.is_gap) {
                            seam = "+" + seam;  // gap
                        } else {
                            seam = "-" + seam;  // overlap
                        }
                    }
                    seam_parts.push_back(seam);
                }
            }
            
            // join seam parts with ":"
            if (!seam_parts.empty()) {
                event.read_seams = "";
                for (size_t i = 0; i < seam_parts.size(); ++i) {
                    if (i > 0) {
                        event.read_seams += ":";
                    }
                    event.read_seams += seam_parts[i];
                }
                resolved_events++;
            }
        }
    }
    
    if (events_with_seams == 0) {
        cout << "no read intervals to resolve for library " << lib_id << endl;
    } else {
        cout << "resolved seams for " << resolved_events << "/" << events_with_seams << " events" << endl;
    }
}

string RearrangeReadEvent::extract_sequence_interval(const string& sequence, const SeamInterval& interval) const
{
    // for overlaps, swap start and end to extract the correct interval
    uint32_t start = interval.start;
    uint32_t end = interval.end;
    if (!interval.is_gap) {
        std::swap(start, end);
    }
    
    // return empty string if no sequence to extract
    if (start == end) {
        return "";
    }
    
    // validate coordinates with assertions
    massert(start < end, "invalid interval: start (%u) >= end (%u) for read %s", 
           start, end, interval.read_id.c_str());
    massert(end <= sequence.length(), "interval end (%u) exceeds sequence length (%zu) for read %s", 
           end, sequence.length(), interval.read_id.c_str());
    
    // extract substring using corrected start and end
    string result = sequence.substr(start, end - start);
    
    // apply reverse complement if needed
    if (interval.reverse_complement) {
        result = reverse_complement(result);
    }
    
    return result;
}


string RearrangeReadEvent::extract_seam(const string& sequence, uint32_t end_left, uint32_t start_right, bool use_placeholders) const
{
    string result = "";
    if (end_left == start_right) {
        return result;
    } else if (end_left < start_right) {
        // gap between alignments
        if (use_placeholders) {
            result = "+" + string(start_right - end_left, 'N');
        } else {
            result = "+" + safe_substr(sequence, end_left, start_right);
        }
    } else {
        // overlap between alignments (end_left > start_right)
        if (use_placeholders) {
            result = "-" + string(end_left - start_right, 'N');
        } else {
            result = "-" + safe_substr(sequence, start_right, end_left);
        }
    }
    return result;
}

string RearrangeReadEvent::extract_assembly_seam(const Alignment& align1, const Alignment& align2) const
{
    // assembly coordinates
    uint32_t end_left = align1.contig_end;
    uint32_t start_right = align2.contig_start;
    massert(align1.contig_start < align2.contig_end, "align1.contig_start must be less than align2.contig_end");

    int32_t margin = abs(static_cast<int32_t>(end_left) - static_cast<int32_t>(start_right));
    if (margin > max_margin) {
        cerr << "assembly seam exceeds max_margin: " << margin << " > " << max_margin << endl;
        print_alignment(align1, store, "align1");
        print_alignment(align2, store, "align2");
        cerr << "end_left=" << end_left << " start_right=" << start_right << endl;
    }
    massert(margin <= max_margin, "assembly seam exceeds max_margin");
    
    string sequence;
    bool use_placeholders = true;
    
    // get contig sequence based on resolve_seams setting
    if (resolve_seams == ResolveSeams::REFERENCE_ONLY || resolve_seams == ResolveSeams::COMPLETE) {
        massert(assembly_sequences != nullptr, "assembly_sequences must be provided for reference seam extraction");
        string contig_id = store.get_contig_id(align1.contig_index);
        sequence = assembly_sequences->get_sequence(contig_id);
        use_placeholders = false;
    }
    
    return extract_seam(sequence, end_left, start_right, use_placeholders);
}

vector<SeamInterval> RearrangeReadEvent::get_seam_intervals(const Alignment& A, const Alignment& B, const Alignment* X,
                                                           RearrangementType type, bool is_read_reversed, const string& lib_id) const
{
    vector<SeamInterval> intervals;
    
    switch (type) {
        case RearrangementType::LARGE_DELETE:
            // delete: read_seams = A:B
            intervals.push_back(create_seam_interval(A, B, lib_id, is_read_reversed));
            break;
            
        case RearrangementType::LARGE_INSERT:
            // insert: read_seams = A:X:B
            intervals.push_back(create_seam_interval(A, *X, lib_id, is_read_reversed));
            intervals.push_back(create_seam_interval(*X, B, lib_id, is_read_reversed));
            break;
            
        case RearrangementType::LARGE_INVERT:
            // invert: read_seams = A:X:B
            intervals.push_back(create_seam_interval(A, *X, lib_id, is_read_reversed));
            intervals.push_back(create_seam_interval(*X, B, lib_id, is_read_reversed));
            break;
            
        default:
            // unknown type
            break;
    }
    
    return intervals;
}

string RearrangeReadEvent::get_assembly_seams(const Alignment& A, const Alignment& B, const Alignment* X,
                                            RearrangementType type) const
{
    switch (type) {
        case RearrangementType::LARGE_DELETE:
            // delete: assembly_seams = empty
            return "";
            
        case RearrangementType::LARGE_INSERT:
            // insert: assembly_seams = A:B
            return extract_assembly_seam(A, B);
            
        case RearrangementType::LARGE_INVERT:
        {
            // invert: assembly_seams = A:X:B
            string assembly_seam1 = extract_assembly_seam(A, *X);
            string assembly_seam2 = extract_assembly_seam(*X, B);
            return assembly_seam1 + ":" + assembly_seam2;
        }
            
        default:
            // unknown type
            return "";
    }
}

ReadEvent RearrangeReadEvent::create_event(const Alignment& A, const Alignment& B, const Alignment* X,
                                          RearrangementType type, bool is_read_reversed, bool verbose) const
{
    if (verbose) {
        cout << "creating event with alignments:" << endl;
        print_alignment(A, store, "  A");
        print_alignment(B, store, "  B");
        if (X) {
            print_alignment(*X, store, "  X");
        }
        cout << "  type: " << type << ", is_read_reversed: " << (is_read_reversed ? "true" : "false") << endl;
    }
    
    // assert X pointer consistency with event type
    if (type == RearrangementType::LARGE_DELETE) {
        massert(X == nullptr, "for delete events, X must be null");
    } else {
        massert(X != nullptr, "for non-delete events, X must not be null");
        // for non-delete events, validate X contig index is reasonable
        // contig indices should be much smaller than this (typically < 100k)
        massert(X->contig_index < 100000, "X contig index out of reasonable bounds: %zu", X->contig_index);
    }
    
    // validate A and B contig indices as well
    massert(A.contig_index < 100000, "A contig index out of reasonable bounds: %zu", A.contig_index);
    massert(B.contig_index < 100000, "B contig index out of reasonable bounds: %zu", B.contig_index);
    
    ReadEvent event;
    
    // get read ID
    event.read_id = store.get_read_id(A.read_index);
    event.type = type;
    event.read_strand = alignment_to_strand_string(is_read_reversed);
    
    // anchor contig information
    event.contig_id = store.get_contig_id(A.contig_index);
    
    // determine clip positions (breakpoints)
    event.out_clip = min(A.contig_end, B.contig_end);
    event.in_clip = max(A.contig_start, B.contig_start);
    event.read_clip_out = min(A.read_end, B.read_end);
    event.read_clip_in = max(A.read_start, B.read_start);
    
    // span coordinates
    event.span_start = A.contig_start;
    event.span_end = B.contig_end;
    event.read_span_start = A.read_start;
    event.read_span_end = B.read_end;
    
    // element information
    if (X) {
        event.element_contig = store.get_contig_id(X->contig_index);
        event.element_strand = alignment_to_strand_string(X->is_reverse);
        event.element_start = X->contig_start;
        event.element_end = X->contig_end;
        
    } else {
        // deletion - no element
        event.element_contig = "";
        event.element_strand = "";
        event.element_start = 0;
        event.element_end = 0;
        
    }
    
    // create read intervals and assembly seams
    event.read_seam_intervals = get_seam_intervals(A, B, X, type, is_read_reversed, current_lib_id);
    event.assembly_seams = get_assembly_seams(A, B, X, type);
    
    // store alignment copies for verification
    event.alignment_A = A;  // copy alignment A
    event.alignment_B = B;  // copy alignment B
    if (X) {
        event.has_alignment_X = true;
        event.alignment_X = *X;  // copy dereferenced alignment X
    } else {
        event.has_alignment_X = false;
        // alignment_X doesn't need to be set when has_alignment_X is false
    }
    
    return event;
}
