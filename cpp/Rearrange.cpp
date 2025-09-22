#include "Rearrange.h"
#include "VerifyRearrange.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

Rearrange::Rearrange(const AlignmentStore& store,
                     VerifyRearrange* verifier,
                     int max_gap,
                     int min_element_length,
                     int min_anchor_length,
                     double max_anchor_mutations_percent,
                     double max_element_mutation_percent)
    : store(store)
    , max_gap(max_gap)
    , min_element_length(min_element_length)
    , min_anchor_length(min_anchor_length)
    , max_anchor_mutations_percent(max_anchor_mutations_percent)
    , max_element_mutation_percent(max_element_mutation_percent)
    , verifier(verifier)
{
}

void Rearrange::get_events()
{
    cout << "detecting rearrangement events..." << endl;
    
    const auto& reads = store.get_reads();
    const auto& alignments = store.get_alignments();
    
    // determine which reads to process
    std::set<uint32_t> reads_to_process;
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
    
    size_t tested_reads = 0;
    size_t reads_with_events = 0;
    size_t found_events = 0;
    
    // iterate through relevant reads only
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
        if (process_read(alignment_indices, alignments, tested_reads, found_events)) {
            reads_with_events++;
        }
    }
    
    cout << "tested " << tested_reads << " read alignment combinations" << endl;
    cout << "found " << found_events << " events in " << reads_with_events << " reads" << endl;
    
    // print rejection summary
    error_tracker.print_summary();
}

bool Rearrange::process_read(const std::vector<size_t>& alignment_indices, 
                            const std::vector<Alignment>& alignments,
                            size_t& tested_reads, size_t& found_events)
{
    bool found_event_in_read = false;
    
    // pre-processing: calculate valid anchors for all alignments
    vector<bool> is_valid_anchor_vec(alignment_indices.size(), false);
    for (size_t i = 0; i < alignment_indices.size(); ++i) {
        const Alignment& aln = alignments[alignment_indices[i]];
        is_valid_anchor_vec[i] = is_valid_anchor(aln);
    }
    
    // DEBUG: hardcoded read name for debug output
    const std::string debug_read = "m84221_250708_233611_s3/136906147/ccs";
    
    // DEBUG: filter for specific read and dump all alignments
    if (!alignment_indices.empty()) {
        const Alignment& first_aln = alignments[alignment_indices[0]];
        string read_id = store.get_reads()[first_aln.read_index].id;
        
        // Only debug specific read
        if (read_id == debug_read) {
            cout << "DEBUG: processing target read " << read_id << endl;
            cout << "DEBUG: dumping all " << alignment_indices.size() << " alignments in order:" << endl;
            
            for (size_t i = 0; i < alignment_indices.size(); ++i) {
                const Alignment& aln = alignments[alignment_indices[i]];
                bool is_valid = is_valid_anchor_vec[i];
                cout << "DEBUG: alignment[" << i << "] idx=" << alignment_indices[i] 
                     << " read:" << aln.read_start << "-" << aln.read_end 
                     << " contig:" << aln.contig_start << "-" << aln.contig_end
                     << " valid_anchor:" << (is_valid ? "YES" : "NO") << endl;
            }
            
            size_t valid_anchor_count = std::count(is_valid_anchor_vec.begin(), is_valid_anchor_vec.end(), true);
            cout << "DEBUG: valid anchors: " << valid_anchor_count << ", total alignments: " << alignment_indices.size() << endl;
        }
    }
    
    // determine if this is the debug read
    bool is_verbose = false;
    if (!alignment_indices.empty()) {
        const Alignment& first_aln = alignments[alignment_indices[0]];
        string read_id = store.get_reads()[first_aln.read_index].id;
        is_verbose = (read_id == debug_read);
    }
    
    // process trio events
    if (process_read_trios(alignment_indices, alignments, is_valid_anchor_vec, tested_reads, found_events, is_verbose)) {
        found_event_in_read = true;
    }
    
    // process pair events
    if (process_read_pairs(alignment_indices, alignments, is_valid_anchor_vec, tested_reads, found_events, is_verbose)) {
        found_event_in_read = true;
    }
    
    return found_event_in_read;
}

bool Rearrange::process_read_trios(const std::vector<size_t>& alignment_indices,
                                  const std::vector<Alignment>& alignments,
                                  const std::vector<bool>& is_valid_anchor_vec,
                                  size_t& tested_reads, size_t& found_events, bool is_verbose)
{
    bool found_event_in_read = false;
    
    // test trio events
    for (size_t i = 0; i < alignment_indices.size(); ++i) {
        // anchor1: current alignment, must be valid anchor
        if (!is_valid_anchor_vec[i]) {
            if (is_verbose) cout << "trio i=" << i << " skip: anchor1 invalid" << endl;
            continue;
        }
        const Alignment& anchor1 = alignments[alignment_indices[i]];
        
        // find anchor2: first valid anchor after anchor1, skipping the very next alignment
        size_t j = i + 2; // skip i+1
        while (j < alignment_indices.size() && !is_valid_anchor_vec[j]) {
            j++;
        }
        if (j >= alignment_indices.size()) {
            if (is_verbose) cout << "trio i=" << i << " skip: no anchor2" << endl;
            continue;
        }
        const Alignment& anchor2 = alignments[alignment_indices[j]];
                 
        // verify that anchors are on same strand and contig
        if (anchor1.is_reverse != anchor2.is_reverse || anchor1.contig_index != anchor2.contig_index) {
            if (is_verbose) cout << "trio i=" << i << " skip: strand/contig mismatch" << endl;
            continue;
        }
        
        // find element: best alignment in interval between anchor1 and anchor2 (with 1bp margin)
        uint32_t interval_start = anchor1.read_end;
        uint32_t interval_end = anchor2.read_start;
        
        const Alignment* element = store.get_alignment_in_read(anchor1.read_index, interval_start, interval_end);
        if (!element || !is_valid_element(*element)) {
            if (is_verbose) cout << "trio i=" << i << " skip: no element" << endl;
            continue;
        }
        
        // print alignment indices for verbose mode
        if (is_verbose) {
            size_t element_global_idx = element - &alignments[0];
            size_t element_i = alignment_indices.size(); // default to invalid index
            for (size_t k = 0; k < alignment_indices.size(); ++k) {
                if (alignment_indices[k] == element_global_idx) {
                    element_i = k;
                    break;
                }
            }
            cout << "trio i=" << i << " j=" << j << " element_i=" << element_i << endl;
        }
        
        // test gaps (allowing max_gap overlap/gap with anchor1 and anchor2)
        int32_t delta1 = static_cast<int32_t>(element->read_start) - static_cast<int32_t>(anchor1.read_end);
        int32_t delta2 = static_cast<int32_t>(anchor2.read_start) - static_cast<int32_t>(element->read_end);
        
        if (abs(delta1) > max_gap || abs(delta2) > max_gap) {
            if (is_verbose) cout << "trio i=" << i << " skip: gap too large" << endl;
            continue;
        }
        
        // define A, B, X (after all other validation)
        bool is_read_reversed = anchor1.is_reverse;  // determine flag from anchor1/anchor2 strand
        const Alignment* A = is_read_reversed ? &anchor2 : &anchor1;
        const Alignment* B = is_read_reversed ? &anchor1 : &anchor2;
        const Alignment* X = element;
        
        // test A B order
        if (A->contig_start >= B->contig_start) {
            if (is_verbose) cout << "trio i=" << i << " skip: wrong A B order" << endl;
            continue;
        }
        
        tested_reads++;
        ReadEvent event;
        
        EventTestResult result = test_trio_event(*A, *B, *X, is_read_reversed, event);
        
        if (result != EventTestResult::FOUND_INSERTION && 
            result != EventTestResult::FOUND_INVERSION) {
            if (is_verbose) cout << "trio i=" << i << " skip: event test failed" << endl;
            continue;
        }
            
        // verify event if requested
        if (!verifier || verifier->verify_event(event, *A, *B, X, store)) {
            read_events.push_back(event);
            found_event_in_read = true;
            found_events++;
        }
    }
    
    return found_event_in_read;
}

bool Rearrange::process_read_pairs(const std::vector<size_t>& alignment_indices,
                                  const std::vector<Alignment>& alignments,
                                  const std::vector<bool>& is_valid_anchor_vec,
                                  size_t& tested_reads, size_t& found_events, bool is_verbose)
{
    bool found_event_in_read = false;
    
    // test pair events
    for (size_t i = 0; i + 1 < alignment_indices.size(); ++i) {
        // both alignments must be valid anchors
        if (!is_valid_anchor_vec[i] || !is_valid_anchor_vec[i + 1]) {
            if (is_verbose) cout << "pair i=" << i << " skip: anchors invalid" << endl;
            continue;
        }
        
        const Alignment& anchor1 = alignments[alignment_indices[i]];
        const Alignment& anchor2 = alignments[alignment_indices[i + 1]];
        
        // test strand
        if (anchor1.is_reverse != anchor2.is_reverse || anchor1.contig_index != anchor2.contig_index) {
            if (is_verbose) cout << "pair i=" << i << " skip: strand/contig mismatch" << endl;
            continue;
        }
        
        // print alignment indices for verbose mode
        if (is_verbose) cout << "pair i=" << i << " j=" << (i + 1) << endl;
        
        // test gaps (allowing max_gap overlap/gap between anchors)
        int32_t delta = static_cast<int32_t>(anchor2.read_start) - static_cast<int32_t>(anchor1.read_end);
        
        if (abs(delta) > max_gap) {
            if (is_verbose) cout << "pair i=" << i << " skip: gap too large" << endl;
            continue;
        }
        
        // define A, B (after all other validation)
        bool is_read_reversed = anchor1.is_reverse;  // determine flag from anchor1/anchor2 strand
        const Alignment* A = is_read_reversed ? &anchor2 : &anchor1;
        const Alignment* B = is_read_reversed ? &anchor1 : &anchor2;
        
        // test A B order
        if (A->contig_start >= B->contig_start) {
            if (is_verbose) cout << "pair i=" << i << " skip: wrong A B order" << endl;
            continue;
        }
        
        tested_reads++;
        ReadEvent event;
        
        EventTestResult result = test_pair_event(*A, *B, is_read_reversed, event);
        
        if (result != EventTestResult::FOUND_DELETION) {
            if (is_verbose) cout << "pair i=" << i << " skip: event test failed" << endl;
            continue;
        }
        
        // verify event if requested
        if (!verifier || verifier->verify_event(event, *A, *B, nullptr, store)) {
            read_events.push_back(event);
            found_event_in_read = true;
            found_events++;
        }
    }
    
    return found_event_in_read;
}

void Rearrange::aggregate_events(const std::map<std::string, std::string>& key_to_event_id)
{
    map<string, AggregatedEvent> event_map;
    
    for (const ReadEvent& read_event : read_events) {
        AggregatedEvent agg_event;
        agg_event.type = read_event.type;
        agg_event.contig_id = read_event.contig_id;
        agg_event.strand = read_event.strand;
        agg_event.out_clip = read_event.out_clip;
        agg_event.in_clip = read_event.in_clip;
        agg_event.element_contig = read_event.element_contig;
        agg_event.element_strand = read_event.element_strand;
        agg_event.element_start = read_event.element_start;
        agg_event.element_end = read_event.element_end;
        agg_event.left_shim = read_event.left_shim;
        agg_event.right_shim = read_event.right_shim;
        agg_event.middle_shim = read_event.middle_shim;
        
        string key = agg_event.create_key();
        
        if (event_map.find(key) != event_map.end()) {
            event_map[key].read_count++;
        } else {
            agg_event.read_count = 1;
            event_map[key] = agg_event;
        }
    }
    
    // convert map to vector and assign event IDs from manager
    aggregated_events.clear();
    
    for (const auto& pair : event_map) {
        AggregatedEvent event = pair.second;
        
        // get event ID from manager's key mapping
        auto it = key_to_event_id.find(pair.first);
        if (it != key_to_event_id.end()) {
            event.event_id = it->second;
        } else {
            // this shouldn't happen if manager is working correctly
            cerr << "warning: event key not found in manager mapping: " << pair.first << endl;
            event.event_id = "unknown";
        }
        
        aggregated_events.push_back(event);
    }
    
    // sort by event ID for consistent ordering
    sort(aggregated_events.begin(), aggregated_events.end(),
         [](const AggregatedEvent& a, const AggregatedEvent& b) {
             return a.event_id < b.event_id;
         });
}

void Rearrange::write_to_csv(const string& ofn_prefix)
{
    // write rejection details
    error_tracker.write_to_file(ofn_prefix);
    
    // write read events
    string read_events_file = ofn_prefix + "_read_events.tsv";
    cout << "writing read events to " << read_events_file << endl;
    
    ofstream ofs_read(read_events_file);
    if (!ofs_read.is_open()) {
        cerr << "error: could not open file " << read_events_file << endl;
        return;
    }
    
    ofs_read << "read_id\ttype\tcontig_id\tstrand\tout_clip\tin_clip\tread_clip_out\tread_clip_in\telement_contig\telement_strand\telement_start\telement_end\tleft_shim\tright_shim\tmiddle_shim\n";
    
    for (const ReadEvent& event : read_events) {
        ofs_read << event.read_id << "\t"
                 << event.type << "\t"
                 << event.contig_id << "\t"
                 << event.strand << "\t"
                 << event.out_clip << "\t"
                 << event.in_clip << "\t"
                 << event.read_clip_out << "\t"
                 << event.read_clip_in << "\t"
                 << event.element_contig << "\t"
                 << event.element_strand << "\t"
                 << event.element_start << "\t"
                 << event.element_end << "\t"
                 << event.left_shim << "\t"
                 << event.right_shim << "\t"
                 << event.middle_shim << "\n";
    }
    ofs_read.close();
    cout << "wrote " << read_events.size() << " read events to " << read_events_file << endl;
    
    // write aggregated events
    string agg_events_file = ofn_prefix + "_aggregated_events.tsv";
    cout << "writing aggregated events to " << agg_events_file << endl;
    
    ofstream ofs_agg(agg_events_file);
    if (!ofs_agg.is_open()) {
        cerr << "error: could not open file " << agg_events_file << endl;
        return;
    }
    
    ofs_agg << "event_id\ttype\tcontig_id\tstrand\tout_clip\tin_clip\telement_contig\telement_strand\telement_start\telement_end\tread_count\tread_coverage\n";
    
    for (const AggregatedEvent& event : aggregated_events) {
        ofs_agg << event.event_id << "\t"
                << event.type << "\t"
                << event.contig_id << "\t"
                << event.strand << "\t"
                << event.out_clip << "\t"
                << event.in_clip << "\t"
                << event.element_contig << "\t"
                << event.element_strand << "\t"
                << event.element_start << "\t"
                << event.element_end << "\t"
                << event.read_count << "\t"
                << event.read_coverage << "\n";
    }
    ofs_agg.close();
    cout << "wrote " << aggregated_events.size() << " aggregated events to " << agg_events_file << endl;
    
    // write read support table
    string read_support_file = ofn_prefix + "_read_support.tsv";
    cout << "writing read support to " << read_support_file << endl;
    
    ofstream ofs_support(read_support_file);
    if (!ofs_support.is_open()) {
        cerr << "error: could not open file " << read_support_file << endl;
        return;
    }
    
    ofs_support << "event_id\tread_id\tcontig\tclip_out\tclip_in\tsupport_type\tdistance\t"
                << "alignment_A_start\talignment_A_end\talignment_A_strand\t"
                << "read_alignment_A_start\tread_alignment_A_end\t"
                << "alignment_B_start\talignment_B_end\talignment_B_strand\t"
                << "read_alignment_B_start\tread_alignment_B_end\tcontains_event\n";
    
    for (const ReadSupport& support : read_supports) {
        ofs_support << support.event_id << "\t"
                    << support.read_id << "\t"
                    << support.contig << "\t"
                    << support.clip_out << "\t"
                    << support.clip_in << "\t"
                    << support.support_type << "\t"
                    << support.distance << "\t"
                    << support.alignment_A_start << "\t"
                    << support.alignment_A_end << "\t"
                    << support.alignment_A_strand << "\t"
                    << support.read_alignment_A_start << "\t"
                    << support.read_alignment_A_end << "\t"
                    << support.alignment_B_start << "\t"
                    << support.alignment_B_end << "\t"
                    << support.alignment_B_strand << "\t"
                    << support.read_alignment_B_start << "\t"
                    << support.read_alignment_B_end << "\t"
                    << (support.contains_event ? "true" : "false") << "\n";
    }
    ofs_support.close();
    cout << "wrote " << read_supports.size() << " read support entries to " << read_support_file << endl;
}

string Rearrange::alignment_to_strand_string(bool is_reverse) const
{
    return is_reverse ? "-" : "+";
}

// new event testing functions implementation
EventTestResult Rearrange::test_trio_event(const Alignment& A, const Alignment& B, const Alignment& X,
                                           bool is_read_reversed, ReadEvent& event)
{
    // assert A is left of B on contig
    if (A.contig_start >= B.contig_start) {
        print_event_debug("ASSERTION FAILURE in test_trio_event", A, X, B, &A, &B, &X, is_read_reversed);
        cerr << "error: assertion failed - A must be left of B" << endl;
        exit(1);
    }
    
    // assert alignments are properly ordered
    if (A.contig_start >= A.contig_end || B.contig_start >= B.contig_end || X.contig_start >= X.contig_end) {
        print_event_debug("ASSERTION FAILURE in test_trio_event - malformed coordinates", A, X, B, &A, &B, &X, is_read_reversed);
        cerr << "error: assertion failed - alignment coordinates malformed" << endl;
        exit(1);
    }
    
    // check for alignment overlaps
    if (alignments_overlap(A, B)) {
        error_tracker.record_rejection(store.get_read_id(A.read_index), 
                                     EventTestResult::REJECT_ALIGNMENT_OVERLAP, A, B, store, &X);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    if (alignments_overlap(A, X)) {
        error_tracker.record_rejection(store.get_read_id(A.read_index), 
                                     EventTestResult::REJECT_ALIGNMENT_OVERLAP, A, B, store, &X);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    if (alignments_overlap(X, B)) {
        error_tracker.record_rejection(store.get_read_id(A.read_index), 
                                     EventTestResult::REJECT_ALIGNMENT_OVERLAP, A, B, store, &X);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    // classify element position and strand
    if (X.contig_end <= A.contig_start || X.contig_start >= B.contig_end || X.contig_index != A.contig_index) {
        // element is outside span → insertion
        event = create_event(A, B, &X, RearrangementType::LARGE_INSERT, is_read_reversed);
        return EventTestResult::FOUND_INSERTION;
    }
    
    // element is inside span - calculate max distance between X and A, and X and B
    uint32_t distance_X_A = X.contig_start - A.contig_end;
    uint32_t distance_X_B = B.contig_start - X.contig_end;
    uint32_t max_distance = std::max(distance_X_A, distance_X_B);
    
    if (max_distance > static_cast<uint32_t>(max_gap)) {
        error_tracker.record_rejection(store.get_read_id(A.read_index), 
                                     EventTestResult::REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE, A, B, store, &X);
        return EventTestResult::REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE;
    }
    
    // check strand orientation
    if (X.is_reverse == is_read_reversed) {
        // X is on same strand as read orientation → insertion
        event = create_event(A, B, &X, RearrangementType::LARGE_INSERT, is_read_reversed);
        return EventTestResult::FOUND_INSERTION;
    } else {
        // X is on different strand from read orientation → inversion
        event = create_event(A, B, &X, RearrangementType::LARGE_INVERT, is_read_reversed);
        return EventTestResult::FOUND_INVERSION;
    }
}

EventTestResult Rearrange::test_pair_event(const Alignment& A, const Alignment& B, 
                                           bool is_read_reversed, ReadEvent& event)
{
    // A, B are already assigned by caller
    
    // assert A is left of B on contig
    if (A.contig_start >= B.contig_start) {
        print_pair_debug("ASSERTION FAILURE in test_pair_event", A, B, &A, &B, is_read_reversed);
        cerr << "error: assertion failed - A must be left of B" << endl;
        exit(1);
    }
    
    // assert alignments are properly ordered
    if (A.contig_start >= A.contig_end || B.contig_start >= B.contig_end) {
        print_pair_debug("ASSERTION FAILURE in test_pair_event - malformed coordinates", A, B, &A, &B, is_read_reversed);
        cerr << "error: assertion failed - alignment coordinates malformed" << endl;
        exit(1);
    }
    
    // check for alignment overlap
    if (alignments_overlap(A, B)) {
        error_tracker.record_rejection(store.get_read_id(A.read_index), 
                                     EventTestResult::REJECT_ALIGNMENT_OVERLAP, A, B, store);
        return EventTestResult::REJECT_ALIGNMENT_OVERLAP;
    }
    
    // check gap in read coords between A and B
    uint32_t gap_A_B = (B.read_start > A.read_end) ? (B.read_start - A.read_end) : 0;
    if (gap_A_B > static_cast<uint32_t>(max_gap)) {
        error_tracker.record_rejection(store.get_read_id(A.read_index), 
                                     EventTestResult::REJECT_READ_GAP_TOO_LARGE, A, B, store);
        return EventTestResult::REJECT_READ_GAP_TOO_LARGE;
    }
    
    // check if inferred element X is long enough
    uint32_t element_length = B.contig_start - A.contig_end;
    if (static_cast<int>(element_length) < min_element_length) {
        error_tracker.record_rejection(store.get_read_id(A.read_index), 
                                     EventTestResult::REJECT_ELEMENT_TOO_SHORT, A, B, store);
        return EventTestResult::REJECT_ELEMENT_TOO_SHORT;
    }
    
    // deletion event
    event = create_event(A, B, nullptr, RearrangementType::LARGE_DELETE, is_read_reversed);
    return EventTestResult::FOUND_DELETION;
}

// helper function implementations
std::string Rearrange::extract_read_shim(const std::string& read_id, uint32_t start, uint32_t end, bool is_read_reversed) const
{
    if (start >= end) {
        return "";  // empty shim
    }
    
    // if no read sequences loaded, return N's as placeholder
    if (read_sequences.empty()) {
        std::string shim(end - start, 'N');
        return shim;
    }
    
    // look up read sequence
    auto read_it = read_sequences.find(read_id);
    if (read_it == read_sequences.end()) {
        // read not found, return N's as placeholder
        std::string shim(end - start, 'N');
        return shim;
    }
    
    const std::string& read_seq = read_it->second;
    
    // check bounds
    if (end > read_seq.length()) {
        // coordinates out of bounds, return N's as placeholder
        std::string shim(end - start, 'N');
        return shim;
    }
    
    // extract shim sequence
    std::string shim = read_seq.substr(start, end - start);
    
    // reverse complement if read is reversed
    if (is_read_reversed) {
        shim = reverse_complement(shim);
    }
    
    return shim;
}

bool Rearrange::alignments_overlap(const Alignment& A, const Alignment& B) const
{
    return !(A.contig_end <= B.contig_start || B.contig_end <= A.contig_start);
}

bool Rearrange::is_valid_element(const Alignment& aln) const
{
    // check minimum length
    uint32_t alignment_length = aln.read_end - aln.read_start;
    if (static_cast<int>(alignment_length) < min_element_length) {
        return false;
    }
    
    // check mutation percentage
    double mutation_percent = static_cast<double>(aln.get_mutation_count()) / static_cast<double>(alignment_length) * 100.0;
    if (mutation_percent > max_element_mutation_percent * 100.0) {
        return false;
    }
    
    return true;
}

bool Rearrange::is_valid_anchor(const Alignment& aln) const
{
    // check minimum length
    uint32_t alignment_length = aln.read_end - aln.read_start;
    if (static_cast<int>(alignment_length) < min_anchor_length) {
        return false;
    }
    
    // check mutation percentage
    double mutation_percent = static_cast<double>(aln.get_mutation_count()) / static_cast<double>(alignment_length) * 100.0;
    if (mutation_percent > max_anchor_mutations_percent * 100.0) {
        return false;
    }
    
    // check alignment overlap
    if (store.get_alignment_overlap(aln) > static_cast<uint32_t>(max_gap)) {
        return false;
    }
    
    // check if alignment is part of relevant alignments
    if (relevant_alignments.find(&aln) == relevant_alignments.end()) {
        return false;
    }
    
    return true;
}

ReadEvent Rearrange::create_event(const Alignment& A, const Alignment& B, const Alignment* X, 
                                  RearrangementType type, bool is_read_reversed) const
{
    ReadEvent event;
    event.read_id = store.get_read_id(A.read_index);
    event.type = type;
    event.contig_id = store.get_contig_id(A.contig_index);
    
    // 1. The end of A is always the clip_out and the beginning of B is always the clip_in
    event.strand = is_read_reversed ? "-" : "+";
    event.out_clip = A.contig_end;
    event.in_clip = B.contig_start;
    event.read_clip_out = A.read_end;
    event.read_clip_in = B.read_start;
    
    // 2. If X is not null then add it to the event
    if (X) {
        event.element_contig = store.get_contig_id(X->contig_index);
        event.element_strand = X->is_reverse ? "-" : "+";
        event.element_start = X->contig_start;
        event.element_end = X->contig_end;
        
        // 3. Assert element_start < element_end (do not swap)
        massert(event.element_start < event.element_end, 
                "element_start (%u) must be < element_end (%u)", 
                event.element_start, event.element_end);
    }
    
    // 4. Extract shims: middle_shim for deletions, left and right for insert/invert
    if (type == RearrangementType::LARGE_DELETE) {
        // for deletions: middle_shim between A and B, left and right are empty
        event.left_shim = "";
        event.right_shim = "";
        event.middle_shim = extract_read_shim(event.read_id, event.read_clip_out, event.read_clip_in, is_read_reversed);
    } else {
        massert(X != nullptr, "X must not be null for insertions and inversions");
        // for insertions and inversions: left_shim (A to X) and right_shim (X to B), middle is empty
        event.left_shim = extract_read_shim(event.read_id, A.read_end, X->read_start, is_read_reversed);
        event.right_shim = extract_read_shim(event.read_id, X->read_end, B.read_start, is_read_reversed);
        event.middle_shim = "";
    }
    
    return event;
}

void Rearrange::get_coverage()
{
    cout << "calculating read coverage for " << aggregated_events.size() << " events..." << endl;
    
    // for each event, calculate coverage by counting unique reads in breakpoint intervals
    for (AggregatedEvent& event : aggregated_events) {
        
        // skip if contig doesn't exist in this store
        if (!store.has_contig_id(event.contig_id)) {
            event.read_coverage = 0;
            continue;
        }
        
        // create two intervals: [clip_out-1, clip_out] and [clip_in, clip_in+1]
        uint32_t clip_out_0based = event.out_clip - 1; // convert to 0-based
        uint32_t clip_in_0based = event.in_clip - 1;   // convert to 0-based
        
        Interval interval1(event.contig_id, 
                          clip_out_0based > 0 ? clip_out_0based - 1 : 0, 
                          clip_out_0based);
        Interval interval2(event.contig_id, clip_in_0based + 1, clip_in_0based + 2);
        
        // get alignments from both intervals
        auto alignments1 = store.get_alignments_in_interval(interval1);
        auto alignments2 = store.get_alignments_in_interval(interval2);
        
        // collect unique read IDs that pass mutation filter
        set<string> reads1;
        set<string> reads2;
        
        // process alignments from first interval
        for (const auto& alignment_ref : alignments1) {
            const Alignment& aln = alignment_ref.get();
            
            // apply mutation filter
            double mutation_percent = static_cast<double>(aln.get_mutation_count()) / 
                static_cast<double>(aln.read_end - aln.read_start);
            if (mutation_percent <= max_anchor_mutations_percent) {
                string read_id = store.get_read_id(aln.read_index);
                reads1.insert(read_id);
            }
        }
        
        // process alignments from second interval
        for (const auto& alignment_ref : alignments2) {
            const Alignment& aln = alignment_ref.get();
            
            // apply mutation filter
            double mutation_percent = static_cast<double>(aln.get_mutation_count()) / 
                static_cast<double>(aln.read_end - aln.read_start);
            if (mutation_percent <= max_anchor_mutations_percent) {
                string read_id = store.get_read_id(aln.read_index);
                reads2.insert(read_id);
            }
        }
        
        // find intersection of read sets
        set<string> intersection;
        set_intersection(reads1.begin(), reads1.end(),
                        reads2.begin(), reads2.end(),
                        inserter(intersection, intersection.begin()));
        
        // set coverage as size of intersection
        event.read_coverage = static_cast<uint32_t>(intersection.size());
    }
    
    cout << "calculated coverage for " << aggregated_events.size() << " events" << endl;
}

std::set<std::string> Rearrange::get_event_keys() const
{
    std::set<std::string> keys;
    for (const ReadEvent& event : read_events) {
        AggregatedEvent temp;
        temp.type = event.type;
        temp.contig_id = event.contig_id;
        temp.strand = event.strand;
        temp.out_clip = event.out_clip;
        temp.in_clip = event.in_clip;
        temp.element_contig = event.element_contig;
        temp.element_strand = event.element_strand;
        temp.element_start = event.element_start;
        temp.element_end = event.element_end;
        temp.left_shim = event.left_shim;
        temp.right_shim = event.right_shim;
        temp.middle_shim = event.middle_shim;
        keys.insert(temp.create_key());
    }
    return keys;
}

void Rearrange::set_reads_and_alignments(const std::vector<Interval>& intervals)
{
    cout << "filtering reads and alignments based on " << intervals.size() << " intervals..." << endl;
    
    relevant_read_indices.clear();
    relevant_alignments.clear();
    
    // collect all unique alignments across all intervals
    std::set<const Alignment*> processed_alignments_set;
    
    for (const auto& interval : intervals) {
        // skip if contig doesn't exist in this store
        if (!store.has_contig_id(interval.contig)) {
            continue;
        }
        
        // get alignments overlapping this interval
        std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval);
        
        for (const auto& alignment_ref : alignments) {
            const auto& aln = alignment_ref.get();
            processed_alignments_set.insert(&aln);
        }
    }
    
    // store relevant alignments and extract read indices
    relevant_alignments = processed_alignments_set;
    for (const Alignment* aln : relevant_alignments) {
        relevant_read_indices.insert(aln->read_index);
    }
    
    cout << "found " << relevant_alignments.size() << " relevant alignments from " 
         << relevant_read_indices.size() << " reads" << endl;
}

const RearrangeErrors& Rearrange::get_error_tracker() const
{
    return error_tracker;
}

void Rearrange::load_read_sequences_from_file(const std::string& filename)
{
    cout << "loading read sequences from file: " << filename << endl;
    
    std::unordered_set<std::string> empty_set;
    read_sequences.clear();
    
    FileType file_type = get_file_type(filename);
    switch (file_type) {
    case FileType::FASTA:
        read_fasta(filename, empty_set, read_sequences);
        break;
    case FileType::FASTQ:
        read_fastq(filename, empty_set, read_sequences);
        break;
    default:
        cerr << "error: unsupported read file type for " << filename << endl;
        exit(1);
    }
    
    cout << "loaded " << read_sequences.size() << " read sequences" << endl;
}

void Rearrange::load_read_sequences_from_map(const std::unordered_map<std::string, std::string>& sequences)
{
    cout << "loading " << sequences.size() << " read sequences from R" << endl;
    read_sequences = sequences;
}

void Rearrange::clear_read_sequences()
{
    read_sequences.clear();
}

void Rearrange::print_event_debug(const std::string& context, const Alignment& S1, const Alignment& S2, 
                                  const Alignment& S3, const Alignment* A, const Alignment* B, 
                                  const Alignment* X, bool is_read_reversed) const
{
    cerr << "=== DEBUG: " << context << " ===" << endl;
    cerr << "is_read_reversed: " << (is_read_reversed ? "true" : "false") << endl;
    
    // print original S1, S2, S3
    cerr << "Original alignments:" << endl;
    cerr << "  S1: contig_idx=" << S1.contig_index << " contig_pos=" << S1.contig_start << "-" << S1.contig_end 
         << " read_pos=" << S1.read_start << "-" << S1.read_end << " strand=" << (S1.is_reverse ? "-" : "+") << endl;
    cerr << "  S2: contig_idx=" << S2.contig_index << " contig_pos=" << S2.contig_start << "-" << S2.contig_end 
         << " read_pos=" << S2.read_start << "-" << S2.read_end << " strand=" << (S2.is_reverse ? "-" : "+") << endl;
    cerr << "  S3: contig_idx=" << S3.contig_index << " contig_pos=" << S3.contig_start << "-" << S3.contig_end 
         << " read_pos=" << S3.read_start << "-" << S3.read_end << " strand=" << (S3.is_reverse ? "-" : "+") << endl;
    
    // print assigned A, B, X
    cerr << "Assigned roles:" << endl;
    if (A) {
        cerr << "  A: contig_idx=" << A->contig_index << " contig_pos=" << A->contig_start << "-" << A->contig_end 
             << " read_pos=" << A->read_start << "-" << A->read_end << " strand=" << (A->is_reverse ? "-" : "+") << endl;
    }
    if (B) {
        cerr << "  B: contig_idx=" << B->contig_index << " contig_pos=" << B->contig_start << "-" << B->contig_end 
             << " read_pos=" << B->read_start << "-" << B->read_end << " strand=" << (B->is_reverse ? "-" : "+") << endl;
    }
    if (X) {
        cerr << "  X: contig_idx=" << X->contig_index << " contig_pos=" << X->contig_start << "-" << X->contig_end 
             << " read_pos=" << X->read_start << "-" << X->read_end << " strand=" << (X->is_reverse ? "-" : "+") << endl;
    }
    
    // print ordering check
    if (A && B) {
        cerr << "Ordering check:" << endl;
        cerr << "  A.contig_start (" << A->contig_start << ") >= B.contig_start (" << B->contig_start << "): " 
             << (A->contig_start >= B->contig_start ? "TRUE (PROBLEM!)" : "false (ok)") << endl;
        cerr << "  A contig coords are properly ordered: " 
             << (A->contig_start < A->contig_end ? "true" : "FALSE (PROBLEM!)") << endl;
        cerr << "  B contig coords are properly ordered: " 
             << (B->contig_start < B->contig_end ? "true" : "FALSE (PROBLEM!)") << endl;
    }
    
    cerr << "===========================================" << endl;
}

void Rearrange::print_pair_debug(const std::string& context, const Alignment& S1, const Alignment& S2,
                                 const Alignment* A, const Alignment* B, bool is_read_reversed) const
{
    cerr << "=== DEBUG: " << context << " ===" << endl;
    cerr << "is_read_reversed: " << (is_read_reversed ? "true" : "false") << endl;
    
    // print original S1, S2
    cerr << "Original alignments:" << endl;
    cerr << "  S1: contig_idx=" << S1.contig_index << " contig_pos=" << S1.contig_start << "-" << S1.contig_end 
         << " read_pos=" << S1.read_start << "-" << S1.read_end << " strand=" << (S1.is_reverse ? "-" : "+") << endl;
    cerr << "  S2: contig_idx=" << S2.contig_index << " contig_pos=" << S2.contig_start << "-" << S2.contig_end 
         << " read_pos=" << S2.read_start << "-" << S2.read_end << " strand=" << (S2.is_reverse ? "-" : "+") << endl;
    
    // print assigned A, B
    cerr << "Assigned roles:" << endl;
    if (A) {
        cerr << "  A: contig_idx=" << A->contig_index << " contig_pos=" << A->contig_start << "-" << A->contig_end 
             << " read_pos=" << A->read_start << "-" << A->read_end << " strand=" << (A->is_reverse ? "-" : "+") << endl;
    }
    if (B) {
        cerr << "  B: contig_idx=" << B->contig_index << " contig_pos=" << B->contig_start << "-" << B->contig_end 
             << " read_pos=" << B->read_start << "-" << B->read_end << " strand=" << (B->is_reverse ? "-" : "+") << endl;
    }
    
    // print ordering check
    if (A && B) {
        cerr << "Ordering check:" << endl;
        cerr << "  A.contig_start (" << A->contig_start << ") >= B.contig_start (" << B->contig_start << "): " 
             << (A->contig_start >= B->contig_start ? "TRUE (PROBLEM!)" : "false (ok)") << endl;
        cerr << "  A contig coords are properly ordered: " 
             << (A->contig_start < A->contig_end ? "true" : "FALSE (PROBLEM!)") << endl;
        cerr << "  B contig coords are properly ordered: " 
             << (B->contig_start < B->contig_end ? "true" : "FALSE (PROBLEM!)") << endl;
    }
    
    cerr << "===========================================" << endl;
}

string AggregatedEvent::create_key() const
{
    ostringstream oss;
    oss << static_cast<int>(type) << "|"
        << contig_id << "|"
        << out_clip << "|"
        << in_clip << "|"
        << element_contig << "|"
        << element_strand << "|"
        << element_start << "|"
        << element_end << "|"
        << left_shim << "|"
        << right_shim << "|"
        << middle_shim;
    return oss.str();
}

