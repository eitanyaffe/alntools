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
    default:
        return "unknown";
    }
}

RearrangeReadEvent::RearrangeReadEvent(const AlignmentStore& store,
                                     RearrangeVerify* verifier,
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
{
}

void RearrangeReadEvent::detect_events(vector<ReadEvent>& output_events, 
                                      map<string, size_t>& output_rejection_counts)
{
    cout << "detecting rearrangement events..." << endl;
    
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
        
        // process this read
        if (process_read(alignment_indices, alignments, output_events, output_rejection_counts)) {
            reads_with_events++;
        }
        
        tested_reads++;
    }
    
    cout << "tested " << tested_reads << " reads" << endl;
    cout << "found " << output_events.size() << " events in " << reads_with_events << " reads" << endl;
    
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

void RearrangeReadEvent::load_read_sequences_from_file(const string& filename)
{
    // TODO: implement FASTA file reading if needed
    cout << "loading read sequences from file: " << filename << " (not implemented)" << endl;
}

void RearrangeReadEvent::load_read_sequences_from_map(const unordered_map<string, string>& sequences)
{
    read_sequences = sequences;
    cout << "loaded " << read_sequences.size() << " read sequences" << endl;
}

void RearrangeReadEvent::clear_read_sequences()
{
    read_sequences.clear();
}

size_t RearrangeReadEvent::get_total_rejected_count(const map<string, size_t>& rejection_counts) const
{
    size_t total = 0;
    for (const auto& pair : rejection_counts) {
        total += pair.second;
    }
    return total;
}

void RearrangeReadEvent::print_event(const ReadEvent& event) const
{
    cout << "DEBUG EVENT: read_id=" << event.read_id 
         << " type=" << static_cast<int>(event.type)
         << " contig=" << event.contig_id 
         << " out_clip=" << event.out_clip 
         << " in_clip=" << event.in_clip
         << " read_strand=" << event.read_strand;
    
    if (!event.element_contig.empty()) {
        cout << " element_contig=" << event.element_contig
             << " element_strand=" << event.element_strand
             << " element_start=" << event.element_start
             << " element_end=" << event.element_end;
        
        if (!event.left_shim.empty()) {
            cout << " left_shim_len=" << event.left_shim.length();
        }
        if (!event.right_shim.empty()) {
            cout << " right_shim_len=" << event.right_shim.length();
        }
    }
    
    if (!event.middle_shim.empty()) {
        cout << " middle_shim_len=" << event.middle_shim.length();
    }
    
    cout << endl;
}

void RearrangeReadEvent::add_rejection(EventTestResult result, map<string, size_t>& rejection_counts)
{
    string reason = event_test_result_to_string(result);
    rejection_counts[reason]++;
}

bool RearrangeReadEvent::is_valid_anchor(const Alignment& aln) const
{
    // check minimum anchor length
    uint32_t alignment_length = aln.read_end - aln.read_start;
    if (static_cast<int>(alignment_length) < min_anchor_length) {
        return false;
    }
    
    // check mutation rate
    double mutation_percent = static_cast<double>(aln.get_mutation_count()) / 
        static_cast<double>(alignment_length);
    if (mutation_percent > max_anchor_mutations_percent) {
        return false;
    }
    
    // check alignment overlap
    if (store.get_alignment_overlap(aln) > static_cast<uint32_t>(max_margin)) {
        return false;
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
                                     map<string, size_t>& rejection_counts)
{
    bool found_event_in_read = false;
    
    // pre-processing: calculate valid anchors for all alignments
    vector<bool> is_valid_anchor_vec(alignment_indices.size(), false);
    for (size_t i = 0; i < alignment_indices.size(); ++i) {
        const Alignment& aln = alignments[alignment_indices[i]];
        is_valid_anchor_vec[i] = is_valid_anchor(aln);
    }
    
    
    // try trio events first (A-X-B pattern)
    if (process_read_trios(alignment_indices, alignments, is_valid_anchor_vec, 
                          output_events, rejection_counts)) {
        found_event_in_read = true;
    }
    
    // try pair events (A-B pattern with gap)
    if (process_read_pairs(alignment_indices, alignments, is_valid_anchor_vec,
                          output_events, rejection_counts)) {
        found_event_in_read = true;
    }
    
    return found_event_in_read;
}

bool RearrangeReadEvent::process_read_trios(const vector<size_t>& alignment_indices,
                                           const vector<Alignment>& alignments,
                                           const vector<bool>& is_valid_anchor_vec,
                                           vector<ReadEvent>& output_events,
                                           map<string, size_t>& rejection_counts)
{
    bool found_event = false;
    
    // test trio events - based on original sequential logic
    for (size_t i = 0; i < alignment_indices.size(); ++i) {
        // anchor1: current alignment, must be valid anchor
        if (!is_valid_anchor_vec[i]) {
            continue;
        }
        const Alignment& anchor1 = alignments[alignment_indices[i]];
        
        // find anchor2: first valid anchor after anchor1, skipping the very next alignment
        size_t j = i + 2; // skip i+1
        while (j < alignment_indices.size() && !is_valid_anchor_vec[j]) {
            j++;
        }
        if (j >= alignment_indices.size()) {
            continue;
        }
        const Alignment& anchor2 = alignments[alignment_indices[j]];
                 
        // verify that anchors are on same strand and contig
        if (anchor1.is_reverse != anchor2.is_reverse || anchor1.contig_index != anchor2.contig_index) {
            continue;
        }
        
        // find element: best alignment in interval between anchor1 and anchor2
        uint32_t interval_start = anchor1.read_end;
        uint32_t interval_end = anchor2.read_start;
        
        const Alignment* element = store.get_alignment_in_read(anchor1.read_index, interval_start, interval_end);
        if (!element || !is_valid_element(*element)) {
            continue;
        }
        
        // test gaps (allowing max_margin overlap/gap with anchor1 and anchor2)
        int32_t delta1 = static_cast<int32_t>(element->read_start) - static_cast<int32_t>(anchor1.read_end);
        int32_t delta2 = static_cast<int32_t>(anchor2.read_start) - static_cast<int32_t>(element->read_end);
        
        if (abs(delta1) > max_margin || abs(delta2) > max_margin) {
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
            continue;
        }
        
        ReadEvent event;
        EventTestResult result = test_trio_event(*A, *B, X, is_read_reversed, event, rejection_counts);
        
        if (result == EventTestResult::FOUND_EVENT) {
            // verify event if requested
            if (!verifier || verifier->verify_event(event)) {
                output_events.push_back(event);
                found_event = true;
            }
        }
    }
    
    return found_event;
}

bool RearrangeReadEvent::process_read_pairs(const vector<size_t>& alignment_indices,
                                           const vector<Alignment>& alignments,
                                           const vector<bool>& is_valid_anchor_vec,
                                           vector<ReadEvent>& output_events,
                                           map<string, size_t>& rejection_counts)
{
    bool found_event = false;
    
    // test pair events - consecutive anchors only
    for (size_t i = 0; i + 1 < alignment_indices.size(); ++i) {
        // both alignments must be valid anchors
        if (!is_valid_anchor_vec[i] || !is_valid_anchor_vec[i + 1]) {
            continue;
        }
        
        const Alignment& anchor1 = alignments[alignment_indices[i]];
        const Alignment& anchor2 = alignments[alignment_indices[i + 1]];
        
        // test strand and contig
        if (anchor1.is_reverse != anchor2.is_reverse || anchor1.contig_index != anchor2.contig_index) {
            continue;
        }
        
        // test gaps (allowing max_margin overlap/gap between anchors)
        int32_t delta = static_cast<int32_t>(anchor2.read_start) - static_cast<int32_t>(anchor1.read_end);
        
        if (abs(delta) > max_margin) {
            continue;
        }
        
        // define A, B (after all other validation)
        bool is_read_reversed = anchor1.is_reverse;  // determine flag from anchor1/anchor2 strand
        const Alignment* A = is_read_reversed ? &anchor2 : &anchor1;
        const Alignment* B = is_read_reversed ? &anchor1 : &anchor2;
        
        // test A B order
        if (A->contig_start >= B->contig_start) {
            continue;
        }
        
        ReadEvent event;
        EventTestResult result = test_pair_event(*A, *B, is_read_reversed, event, rejection_counts);
        
        if (result == EventTestResult::FOUND_EVENT) {
            // verify event if requested
            if (!verifier || verifier->verify_event(event)) {
                output_events.push_back(event);
                found_event = true;
            }
        }
    }
    
    return found_event;
}

EventTestResult RearrangeReadEvent::test_trio_event(const Alignment& A, const Alignment& B, const Alignment& X,
                                                   bool is_read_reversed, ReadEvent& event,
                                                   map<string, size_t>& rejection_counts)
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
    if (alignments_overlap(A, B)) {
        add_rejection(EventTestResult::REJECT_ALIGNMENT_OVERLAP, rejection_counts);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    if (alignments_overlap(A, X)) {
        add_rejection(EventTestResult::REJECT_ALIGNMENT_OVERLAP, rejection_counts);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    if (alignments_overlap(X, B)) {
        add_rejection(EventTestResult::REJECT_ALIGNMENT_OVERLAP, rejection_counts);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    // classify element position and strand
    if (X.contig_end <= A.contig_end || X.contig_start >= B.contig_start || X.contig_index != A.contig_index) {
        // element is outside span → insertion
        event = create_event(A, B, &X, RearrangementType::LARGE_INSERT, is_read_reversed);
        return EventTestResult::FOUND_EVENT;
    }
    
    // element is inside span - calculate max distance between X and A, and X and B
    uint32_t distance_X_A = X.contig_start - A.contig_end;
    uint32_t distance_X_B = B.contig_start - X.contig_end;
    uint32_t max_distance = max(distance_X_A, distance_X_B);
    
    if (max_distance > static_cast<uint32_t>(max_margin)) {
        add_rejection(EventTestResult::REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE, rejection_counts);
        return EventTestResult::REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE;
    }
    
    // check strand orientation
    if (X.is_reverse) {
        // X is on different strand from read orientation → insertion
        event = create_event(A, B, &X, RearrangementType::LARGE_INVERT, is_read_reversed);
        return EventTestResult::FOUND_EVENT;
    } else {
        // X is on same strand as read orientation → reject as re-insertion
        add_rejection(EventTestResult::REJECT_INSERT_ELEMENT_REINSERTED, rejection_counts);
        return EventTestResult::REJECT_INSERT_ELEMENT_REINSERTED;
    }
}

EventTestResult RearrangeReadEvent::test_pair_event(const Alignment& A, const Alignment& B,
                                                   bool is_read_reversed, ReadEvent& event,
                                                   map<string, size_t>& rejection_counts)
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
    if (alignments_overlap(A, B)) {
        add_rejection(EventTestResult::REJECT_ALIGNMENT_OVERLAP, rejection_counts);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    // check gap in read coords between A and B
    uint32_t gap_A_B = (B.read_start > A.read_end) ? (B.read_start - A.read_end) : 0;
    if (gap_A_B > static_cast<uint32_t>(max_margin)) {
        add_rejection(EventTestResult::REJECT_READ_GAP_TOO_LARGE, rejection_counts);
        return EventTestResult::REJECT_READ_GAP_TOO_LARGE;
    }
    
    // check if inferred element X is long enough
    uint32_t element_length = B.contig_start - A.contig_end;
    if (static_cast<int>(element_length) < min_element_length) {
        add_rejection(EventTestResult::REJECT_ELEMENT_TOO_SHORT, rejection_counts);
        return EventTestResult::REJECT_ELEMENT_TOO_SHORT;
    }
    
    // deletion event
    event = create_event(A, B, nullptr, RearrangementType::LARGE_DELETE, is_read_reversed);
    return EventTestResult::FOUND_EVENT;
}

string RearrangeReadEvent::extract_read_shim(const string& read_id, uint32_t start, uint32_t end, bool is_read_reversed) const
{
    // bounds checking
    if (start >= end) {
        cerr << "warning: invalid shim coordinates for read " << read_id << ": start=" << start << " end=" << end << endl;
        return "";
    }
    
    uint32_t shim_length = end - start;
    const uint32_t MAX_SHIM_LENGTH = 100000; // 100kb limit
    
    if (shim_length > MAX_SHIM_LENGTH) {
        cerr << "warning: shim too long for read " << read_id << ": " << shim_length << " bp, truncating to " << MAX_SHIM_LENGTH << endl;
        end = start + MAX_SHIM_LENGTH;
        shim_length = MAX_SHIM_LENGTH;
    }
    
    auto it = read_sequences.find(read_id);
    if (it == read_sequences.end()) {
        // no sequence available, return N placeholders
        return string(shim_length, 'N');
    }
    
    const string& sequence = it->second;
    if (end > sequence.length()) {
        cerr << "warning: shim coordinates beyond sequence length for read " << read_id << ": end=" << end << " seq_len=" << sequence.length() << endl;
        return string(shim_length, 'N');
    }
    
    string shim = sequence.substr(start, shim_length);
    
    // TODO: handle reverse complement if needed
    if (is_read_reversed) {
        // implement reverse complement
    }
    
    return shim;
}

string RearrangeReadEvent::alignment_to_strand_string(bool is_reverse) const
{
    return is_reverse ? "-" : "+";
}

bool RearrangeReadEvent::alignments_overlap(const Alignment& A, const Alignment& B) const
{
    // check read coordinate overlap
    if (A.read_start < B.read_end && B.read_start < A.read_end) {
        return true;
    }
    
    // check contig coordinate overlap (if on same contig)
    if (A.contig_index == B.contig_index) {
        if (A.contig_start < B.contig_end && B.contig_start < A.contig_end) {
            return true;
        }
    }
    
    return false;
}

ReadEvent RearrangeReadEvent::create_event(const Alignment& A, const Alignment& B, const Alignment* X,
                                          RearrangementType type, bool is_read_reversed) const
{
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
    
    // element information
    if (X) {
        event.element_contig = store.get_contig_id(X->contig_index);
        event.element_strand = alignment_to_strand_string(X->is_reverse);
        event.element_start = X->contig_start;
        event.element_end = X->contig_end;
        
        // extract shims with debugging
        cout << "extracting shims for read " << event.read_id << ": A.read_end=" << A.read_end 
             << " X.read_start=" << X->read_start << " X.read_end=" << X->read_end 
             << " B.read_start=" << B.read_start << endl;
        
        event.left_shim = extract_read_shim(event.read_id, A.read_end, X->read_start, is_read_reversed);
        event.right_shim = extract_read_shim(event.read_id, X->read_end, B.read_start, is_read_reversed);
    } else {
        // deletion - no element
        event.element_contig = "";
        event.element_strand = "";
        event.element_start = 0;
        event.element_end = 0;
        
        // middle shim for deletion with debugging
        cout << "extracting middle shim for read " << event.read_id << ": A.read_end=" << A.read_end 
             << " B.read_start=" << B.read_start << endl;
        
        event.middle_shim = extract_read_shim(event.read_id, A.read_end, B.read_start, is_read_reversed);
    }
    
    return event;
}
