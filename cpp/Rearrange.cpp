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
                     double max_mutations_percent)
    : store(store)
    , max_gap(max_gap)
    , min_element_length(min_element_length)
    , min_anchor_length(min_anchor_length)
    , max_mutations_percent(max_mutations_percent)
    , verifier(verifier)
{
}

void Rearrange::execute()
{
    cout << "executing rearrangement detection..." << endl;
    
    // ensure read-to-alignment index is built
    if (!store.is_read_alignment_index_built()) {
        cerr << "error: read alignment index not built. call init_read_alignment_index() first" << endl;
        return;
    }
    
    // step 1: detect rearrangement events
    get_events();
    
    // step 2: calculate coverage for detected events
    get_coverage();
    
    cout << "rearrangement detection completed" << endl;
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
    
    size_t processed_reads = 0;
    size_t reads_with_events = 0;
    
    // iterate through relevant reads only
    for (uint32_t read_idx : reads_to_process) {
        
        // get alignments for this read (already sorted by read_start)
        auto it = store.get_read_to_alignment_indices().find(static_cast<uint32_t>(read_idx));
        if (it == store.get_read_to_alignment_indices().end()) {
            continue; // no alignments for this read
        }
        
        const vector<size_t>& alignment_indices = it->second;
        
        // skip reads with less than 3 alignments
        if (alignment_indices.size() < 3) {
            continue;
        }
        
        processed_reads++;
        bool found_event_in_read = false;
        
        // check all consecutive triplets A, X, B
        for (size_t i = 0; i + 2 < alignment_indices.size(); ++i) {
            const Alignment& A = alignments[alignment_indices[i]];
            const Alignment& X = alignments[alignment_indices[i + 1]];
            const Alignment& B = alignments[alignment_indices[i + 2]];
            
            // first check if A and B are valid anchors
            if (!is_valid_anchor(A) || !is_valid_anchor(B)) {
                continue;
            }
            
            // check if A and B are compatible (same contig, same strand)
            if (!are_compatible_anchors(A, B)) {
                continue;
            }
            
            // classify event type based on A-B distance
            RearrangementType type = classify_event(A, B, &X);
            
            // validate geometry for this event type
            if (!validate_event_geometry(A, B, &X, type)) {
                continue;
            }
            
            // create normalized event
            ReadEvent event = create_normalized_event(A, B, &X, type);
            
            // verify event if requested
            if (verifier) {
                if (!verifier->verify_event(event, A, B, &X, store)) {
                    cout << "verification failed for event in read " << event.read_id << " - skipping" << endl;
                    continue;
                }
            }
            
            // final filtering: check if both A and B alignments are relevant (interval filtering)
            if (!relevant_alignments.empty()) {
                if (relevant_alignments.find(&A) == relevant_alignments.end() || 
                    relevant_alignments.find(&B) == relevant_alignments.end()) {
                    continue; // skip event if A or B is not in relevant alignments
                }
            }
            
            read_events.push_back(event);
            found_event_in_read = true;
        }
        
        if (found_event_in_read) {
            reads_with_events++;
        }
    }
    
    cout << "processed " << processed_reads << " reads with >=3 alignments" << endl;
    cout << "found " << read_events.size() << " events in " << reads_with_events << " reads" << endl;
    
    // aggregate events
    aggregate_events();
    
    cout << "aggregated into " << aggregated_events.size() << " unique events" << endl;
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
    if (mutation_percent > max_mutations_percent * 100.0) {
        return false;
    }
    
    return true;
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
    if (mutation_percent > max_mutations_percent * 100.0) {
        return false;
    }
    
    return true;
}

bool Rearrange::are_compatible_anchors(const Alignment& A, const Alignment& B) const
{
    // must be on same contig
    if (A.contig_index != B.contig_index) {
        return false;
    }
    
    // must be on same strand
    if (A.is_reverse != B.is_reverse) {
        return false;
    }
    
    return true;
}

RearrangementType Rearrange::classify_event(const Alignment& A, const Alignment& B, const Alignment* X) const
{
    // calculate gap between A and B on read
    uint32_t read_gap;
    if (!A.is_reverse) {
        // forward strand: gap = B.read_start - A.read_end
        read_gap = (B.read_start > A.read_end) ? (B.read_start - A.read_end) : 0;
    } else {
        // reverse strand: gap = A.read_start - B.read_end (A and B swapped in read order)
        read_gap = (A.read_start > B.read_end) ? (A.read_start - B.read_end) : 0;
    }
    
    // classify based on gap size
    if (static_cast<int>(read_gap) <= max_gap) {
        return RearrangementType::LARGE_INSERT;
    } else if (static_cast<int>(read_gap) >= min_element_length) {
        return RearrangementType::LARGE_DELETE;
    } else {
        // intermediate gap size - check if X could be an inversion
        if (X && X->contig_index == A.contig_index && X->is_reverse != A.is_reverse) {
            return RearrangementType::LARGE_INVERT;
        } else {
            return RearrangementType::LARGE_INSERT; // default fallback
        }
    }
}

bool Rearrange::validate_event_geometry(const Alignment& A, const Alignment& B, const Alignment* X, RearrangementType type) const
{
    switch (type) {
    case RearrangementType::LARGE_INSERT:
        // X must be valid element
        if (!X || !is_valid_element(*X)) {
            return false;
        }
        // for insertion, X must not overlap with A-B span on contig (if same contig)
        if (X->contig_index == A.contig_index && elements_overlap_anchor_span(A, B, *X)) {
            return false;
        }
        // gap between A and B must be small
        if (!A.is_reverse) {
            int gap = static_cast<int>(B.read_start) - static_cast<int>(A.read_end);
            if (abs(gap) > max_gap) {
                return false;
            }
        } else {
            int gap = static_cast<int>(A.read_start) - static_cast<int>(B.read_end);
            if (abs(gap) > max_gap) {
                return false;
            }
        }
        break;
        
    case RearrangementType::LARGE_INVERT:
        // X must be valid element
        if (!X || !is_valid_element(*X)) {
            return false;
        }
        // X must be on same contig as A and B
        if (X->contig_index != A.contig_index) {
            return false;
        }
        // X must be on opposite strand from A and B
        if (X->is_reverse == A.is_reverse) {
            return false;
        }
        // X must fit between A and B on contig with max_gap tolerance
        uint32_t out_clip, in_clip;
        if (!A.is_reverse) {
            out_clip = A.contig_end;
            in_clip = B.contig_start;
        } else {
            out_clip = B.contig_end;
            in_clip = A.contig_start;
        }
        if (abs(static_cast<int>(X->contig_start) - static_cast<int>(out_clip)) > max_gap ||
            abs(static_cast<int>(X->contig_end) - static_cast<int>(in_clip)) > max_gap) {
            return false;
        }
        break;
        
    case RearrangementType::LARGE_DELETE:
        // no element required for deletion
        // gap between A and B must be large enough
        if (!A.is_reverse) {
            int gap = static_cast<int>(B.read_start) - static_cast<int>(A.read_end);
            if (gap < min_element_length) {
                return false;
            }
        } else {
            int gap = static_cast<int>(A.read_start) - static_cast<int>(B.read_end);
            if (gap < min_element_length) {
                return false;
            }
        }
        break;
    }
    
    return true;
}

ReadEvent Rearrange::create_normalized_event(const Alignment& A, const Alignment& B, const Alignment* X, RearrangementType type) const
{
    ReadEvent event;
    event.read_id = store.get_read_id(A.read_index);
    event.type = type;
    event.contig_id = store.get_contig_id(A.contig_index);
    
    // normalize coordinates based on strand
    if (!A.is_reverse) {
        // forward strand
        event.strand = "+";
        event.out_clip = A.contig_end;
        event.in_clip = B.contig_start;
        event.read_clip_out = A.read_end;
        event.read_clip_in = B.read_start;
        
        if (X && type != RearrangementType::LARGE_DELETE) {
            event.element_contig = store.get_contig_id(X->contig_index);
            event.element_strand = alignment_to_strand_string(X->is_reverse);
            event.element_start = X->contig_start;
            event.element_end = X->contig_end;
        }
    } else {
        // reverse strand - swap A and B roles, flip element strand
        event.strand = "-";
        event.out_clip = B.contig_end;
        event.in_clip = A.contig_start;
        event.read_clip_out = B.read_end;
        event.read_clip_in = A.read_start;
        
        if (X && type != RearrangementType::LARGE_DELETE) {
            event.element_contig = store.get_contig_id(X->contig_index);
            event.element_strand = alignment_to_strand_string(!X->is_reverse); // flip element strand
            event.element_start = X->contig_start;
            event.element_end = X->contig_end;
        }
    }
    
    // ensure element_start < element_end
    if (event.element_start > event.element_end) {
        swap(event.element_start, event.element_end);
    }
    
    return event;
}

void Rearrange::aggregate_events()
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
        
        string key = agg_event.create_key();
        
        // build event-to-supporting-reads mapping
        event_to_supporting_reads[key].insert(read_event.read_id);
        
        if (event_map.find(key) != event_map.end()) {
            event_map[key].read_count++;
        } else {
            agg_event.read_count = 1;
            event_map[key] = agg_event;
        }
    }
    
    // convert map to vector and assign event IDs
    aggregated_events.clear();
    int event_counter = 1;
    map<string, string> key_to_event_id; // for updating the supporting reads mapping
    
    for (const auto& pair : event_map) {
        AggregatedEvent event = pair.second;
        event.event_id = "ev" + to_string(event_counter++);
        key_to_event_id[pair.first] = event.event_id;
        aggregated_events.push_back(event);
    }
    
    // update event-to-supporting-reads mapping to use event IDs instead of keys
    map<string, set<string>> updated_event_to_supporting_reads;
    for (const auto& pair : event_to_supporting_reads) {
        string event_id = key_to_event_id[pair.first];
        updated_event_to_supporting_reads[event_id] = pair.second;
    }
    event_to_supporting_reads = updated_event_to_supporting_reads;
    
    // sort by read count (descending)
    sort(aggregated_events.begin(), aggregated_events.end(),
         [](const AggregatedEvent& a, const AggregatedEvent& b) {
             return a.read_count > b.read_count;
         });
}

void Rearrange::write_to_csv(const string& ofn_prefix)
{
    // write read events
    string read_events_file = ofn_prefix + "_read_events.tsv";
    cout << "writing read events to " << read_events_file << endl;
    
    ofstream ofs_read(read_events_file);
    if (!ofs_read.is_open()) {
        cerr << "error: could not open file " << read_events_file << endl;
        return;
    }
    
    ofs_read << "read_id\ttype\tcontig_id\tstrand\tout_clip\tin_clip\tread_clip_out\tread_clip_in\telement_contig\telement_strand\telement_start\telement_end\n";
    
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
                 << event.element_end << "\n";
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

bool Rearrange::elements_overlap_anchor_span(const Alignment& A, const Alignment& B, const Alignment& X) const
{
    uint32_t span_start = min(A.contig_start, B.contig_start);
    uint32_t span_end = max(A.contig_end, B.contig_end);
    
    // check if X overlaps with [span_start, span_end]
    return !(X.contig_end <= span_start || X.contig_start >= span_end);
}

void Rearrange::get_coverage()
{
    cout << "calculating read coverage for events..." << endl;
    
    // build contig-to-events mapping for efficient lookup
    build_contig_to_events_mapping();
    
    const auto& reads = store.get_reads();
    const auto& alignments = store.get_alignments();
    
    // track reads that have already supported each event to prevent duplicates
    map<string, set<string>> event_supported_reads;
    
    size_t total_support_found = 0;
    
    // iterate through each read
    for (size_t read_idx = 0; read_idx < reads.size(); ++read_idx) {
        
        // get alignments for this read
        auto it = store.get_read_to_alignment_indices().find(static_cast<uint32_t>(read_idx));
        if (it == store.get_read_to_alignment_indices().end()) {
            continue; // no alignments for this read
        }
        
        const vector<size_t>& alignment_indices = it->second;
        const string& read_id = store.get_read_id(read_idx);
        
        // check all three support scenarios for this read
        for (size_t i = 0; i < alignment_indices.size(); ++i) {
            const Alignment& A = alignments[alignment_indices[i]];
            
            // skip invalid alignments
            if (!is_valid_support_alignment(A)) {
                continue;
            }
            
            string contig_id = store.get_contig_id(A.contig_index);
            
            // scenario 1: single alignment A contains both clip_out and clip_in
            vector<string> candidate_events_A = find_candidate_events(contig_id, 
                A.contig_start - max_gap, A.contig_end + max_gap);
            
            for (const string& event_id : candidate_events_A) {
                // skip if this read already supports this event
                if (event_supported_reads[event_id].count(read_id)) {
                    continue;
                }
                
                // find the event details
                auto event_it = find_if(aggregated_events.begin(), aggregated_events.end(),
                    [&event_id](const AggregatedEvent& e) { return e.event_id == event_id; });
                if (event_it == aggregated_events.end()) continue;
                
                if (check_event_support(event_id, event_it->out_clip, event_it->in_clip, A, nullptr, SupportType::A)) {
                    bool contains_event = event_to_supporting_reads[event_id].count(read_id) > 0;
                    ReadSupport support = create_read_support(event_id, read_id, contig_id, 
                        event_it->out_clip, event_it->in_clip, A, nullptr, SupportType::A, contains_event);
                    read_supports.push_back(support);
                    event_supported_reads[event_id].insert(read_id);
                    total_support_found++;
                    continue; // found support, move to next alignment
                }
            }
            
            // scenario 2: two alignments A-B
            if (i + 1 < alignment_indices.size()) {
                const Alignment& B = alignments[alignment_indices[i + 1]];
                
                if (is_valid_support_alignment(B) && A.contig_index == B.contig_index && 
                    A.is_reverse == B.is_reverse) {
                    
                    uint32_t search_start = min(A.contig_end, B.contig_start) - max_gap;
                    uint32_t search_end = max(A.contig_end, B.contig_start) + max_gap;
                    
                    vector<string> candidate_events_AB = find_candidate_events(contig_id, search_start, search_end);
                    
                    for (const string& event_id : candidate_events_AB) {
                        // skip if this read already supports this event
                        if (event_supported_reads[event_id].count(read_id)) {
                            continue;
                        }
                        
                        auto event_it = find_if(aggregated_events.begin(), aggregated_events.end(),
                            [&event_id](const AggregatedEvent& e) { return e.event_id == event_id; });
                        if (event_it == aggregated_events.end()) continue;
                        
                        if (check_event_support(event_id, event_it->out_clip, event_it->in_clip, A, &B, SupportType::AB)) {
                            bool contains_event = event_to_supporting_reads[event_id].count(read_id) > 0;
                            ReadSupport support = create_read_support(event_id, read_id, contig_id, 
                                event_it->out_clip, event_it->in_clip, A, &B, SupportType::AB, contains_event);
                            read_supports.push_back(support);
                            event_supported_reads[event_id].insert(read_id);
                            total_support_found++;
                            break; // found support, don't check AXB for this A
                        }
                    }
                }
            }
            
            // scenario 3: three alignments A-X-B
            if (i + 2 < alignment_indices.size()) {
                const Alignment& X = alignments[alignment_indices[i + 1]];
                const Alignment& B = alignments[alignment_indices[i + 2]];
                
                if (is_valid_support_alignment(B) && A.contig_index == B.contig_index && 
                    A.is_reverse == B.is_reverse) {
                    
                    // X just needs to pass mutation filter, not length filter
                    double x_mutation_percent = static_cast<double>(X.get_mutation_count()) / 
                        static_cast<double>(X.read_end - X.read_start) * 100.0;
                    if (x_mutation_percent <= max_mutations_percent * 100.0) {
                        
                        uint32_t search_start = min(A.contig_end, B.contig_start) - max_gap;
                        uint32_t search_end = max(A.contig_end, B.contig_start) + max_gap;
                        
                        vector<string> candidate_events_AXB = find_candidate_events(contig_id, search_start, search_end);
                        
                        for (const string& event_id : candidate_events_AXB) {
                            // skip if this read already supports this event
                            if (event_supported_reads[event_id].count(read_id)) {
                                continue;
                            }
                            
                            auto event_it = find_if(aggregated_events.begin(), aggregated_events.end(),
                                [&event_id](const AggregatedEvent& e) { return e.event_id == event_id; });
                            if (event_it == aggregated_events.end()) continue;
                            
                            if (check_event_support(event_id, event_it->out_clip, event_it->in_clip, A, &B, SupportType::AXB)) {
                                bool contains_event = event_to_supporting_reads[event_id].count(read_id) > 0;
                                ReadSupport support = create_read_support(event_id, read_id, contig_id, 
                                    event_it->out_clip, event_it->in_clip, A, &B, SupportType::AXB, contains_event);
                                read_supports.push_back(support);
                                event_supported_reads[event_id].insert(read_id);
                                total_support_found++;
                                break; // found support, don't need to check further
                            }
                        }
                    }
                }
            }
        }
    }
    
    // update read_coverage for each event
    for (AggregatedEvent& event : aggregated_events) {
        event.read_coverage = event_supported_reads[event.event_id].size();
    }
    
    cout << "found " << total_support_found << " read support instances" << endl;
    cout << "calculated coverage for " << aggregated_events.size() << " events" << endl;
}

void Rearrange::build_contig_to_events_mapping()
{
    contig_to_events.clear();
    
    for (const AggregatedEvent& event : aggregated_events) {
        // store event_id and clip coordinates for efficient lookup
        contig_to_events[event.contig_id].push_back(
            make_pair(event.event_id, make_pair(event.out_clip, event.in_clip)));
    }
    
    // sort events by out_clip coordinate for binary search
    for (auto& pair : contig_to_events) {
        sort(pair.second.begin(), pair.second.end(),
             [](const auto& a, const auto& b) {
                 return a.second.first < b.second.first; // sort by out_clip
             });
    }
}

bool Rearrange::is_valid_support_alignment(const Alignment& aln) const
{
    // check minimum length (same as anchor requirements)
    uint32_t alignment_length = aln.read_end - aln.read_start;
    if (static_cast<int>(alignment_length) < min_anchor_length) {
        return false;
    }
    
    // check mutation percentage
    double mutation_percent = static_cast<double>(aln.get_mutation_count()) / static_cast<double>(alignment_length) * 100.0;
    if (mutation_percent > max_mutations_percent * 100.0) {
        return false;
    }
    
    return true;
}

vector<string> Rearrange::find_candidate_events(const string& contig_id, uint32_t start, uint32_t end) const
{
    vector<string> candidates;
    
    auto it = contig_to_events.find(contig_id);
    if (it == contig_to_events.end()) {
        return candidates;
    }
    
    const auto& events = it->second;
    
    // find events whose clip coordinates overlap with [start, end]
    for (const auto& event : events) {
        uint32_t out_clip = event.second.first;
        uint32_t in_clip = event.second.second;
        
        // check if event interval [out_clip, in_clip] overlaps with [start, end]
        if (!(in_clip < start || out_clip > end)) {
            candidates.push_back(event.first);
        }
    }
    
    return candidates;
}

bool Rearrange::check_event_support(const string& /* event_id */, uint32_t clip_out, uint32_t clip_in,
                                   const Alignment& A, const Alignment* B, SupportType support_type) const
{
    switch (support_type) {
    case SupportType::A:
        // single alignment must contain both clip coordinates with margin
        return (A.contig_start <= clip_out + max_gap) && (A.contig_end >= clip_in - max_gap);
        
    case SupportType::AB:
    case SupportType::AXB:
        // both clip coordinates must be between A and B with margin
        if (!B) return false;
        
        uint32_t gap_start = min(A.contig_end, B->contig_start);
        uint32_t gap_end = max(A.contig_end, B->contig_start);
        
        return (clip_out >= gap_start - max_gap) && (clip_out <= A.contig_end + max_gap) &&
               (clip_in <= gap_end + max_gap) && (clip_in >= B->contig_start - max_gap);
    }
    
    return false;
}

ReadSupport Rearrange::create_read_support(const string& event_id, const string& read_id,
                                          const string& contig, uint32_t clip_out, uint32_t clip_in,
                                          const Alignment& A, const Alignment* B, SupportType support_type,
                                          bool contains_event) const
{
    ReadSupport support;
    support.event_id = event_id;
    support.read_id = read_id;
    support.contig = contig;
    support.clip_out = clip_out;
    support.clip_in = clip_in;
    support.support_type = support_type;
    support.contains_event = contains_event;
    
    // alignment A (always present)
    support.alignment_A_start = A.contig_start;
    support.alignment_A_end = A.contig_end;
    support.alignment_A_strand = alignment_to_strand_string(A.is_reverse);
    support.read_alignment_A_start = A.read_start;
    support.read_alignment_A_end = A.read_end;
    
    if (B && (support_type == SupportType::AB || support_type == SupportType::AXB)) {
        // alignment B (for AB and AXB types)
        support.alignment_B_start = B->contig_start;
        support.alignment_B_end = B->contig_end;
        support.alignment_B_strand = alignment_to_strand_string(B->is_reverse);
        support.read_alignment_B_start = B->read_start;
        support.read_alignment_B_end = B->read_end;
        support.distance = abs(static_cast<int>(B->contig_start) - static_cast<int>(A.contig_end));
    } else {
        // no alignment B (for A type)
        support.alignment_B_start = 0;
        support.alignment_B_end = 0;
        support.alignment_B_strand = "none";
        support.read_alignment_B_start = 0;
        support.read_alignment_B_end = 0;
        support.distance = 0;
    }
    
    return support;
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
        keys.insert(temp.create_key());
    }
    return keys;
}

void Rearrange::assign_event_ids(const std::map<std::string, std::string>& key_to_event_id)
{
    // update event-to-supporting-reads mapping to use assigned event IDs
    std::map<std::string, std::set<std::string>> updated_mapping;
    for (const auto& pair : event_to_supporting_reads) {
        auto it = key_to_event_id.find(pair.first);
        if (it != key_to_event_id.end()) {
            updated_mapping[it->second] = pair.second;
        }
    }
    event_to_supporting_reads = updated_mapping;
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

string AggregatedEvent::create_key() const
{
    ostringstream oss;
    oss << static_cast<int>(type) << "|"
        << contig_id << "|"
        << strand << "|"
        << out_clip << "|"
        << in_clip << "|"
        << element_contig << "|"
        << element_strand << "|"
        << element_start << "|"
        << element_end;
    return oss.str();
}

