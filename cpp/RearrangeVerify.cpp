#include "RearrangeVerify.h"
#include "RearrangeReadEvent.h"
#include "utils.h"
#include "alignment_store.h"
#include <fstream>
#include <iostream>

using namespace std;

// Conversion functions for VerifyErrorMode
VerifyErrorMode string_to_verify_error_mode(const string& str)
{
    if (str == "exit_on_error") return VerifyErrorMode::EXIT_ON_ERROR;
    if (str == "exit_if_error") return VerifyErrorMode::EXIT_IF_ERROR;
    if (str == "warning_only") return VerifyErrorMode::WARNING_ONLY;
    
    cerr << "error: unknown verify error mode: " << str << endl;
    exit(1);
}

string verify_error_mode_to_string(VerifyErrorMode mode)
{
    switch (mode) {
        case VerifyErrorMode::EXIT_ON_ERROR: return "exit_on_error";
        case VerifyErrorMode::EXIT_IF_ERROR: return "exit_if_error";
        case VerifyErrorMode::WARNING_ONLY: return "warning_only";
        default: return "unknown";
    }
}

string verify_result_type_to_string(VerifyResultType result_type)
{
    switch (result_type) {
        case VerifyResultType::SUCCESS: return "success";
        case VerifyResultType::FAILED_ASSEMBLY_MUTATION_A: return "failed_assembly_mutation_A";
        case VerifyResultType::FAILED_ASSEMBLY_MUTATION_B: return "failed_assembly_mutation_B";
        case VerifyResultType::FAILED_ASSEMBLY_MUTATION_X: return "failed_assembly_mutation_X";
        case VerifyResultType::FAILED_READ_SEQUENCE_MATCH_A: return "failed_read_sequence_match_A";
        case VerifyResultType::FAILED_READ_SEQUENCE_MATCH_B: return "failed_read_sequence_match_B";
        case VerifyResultType::FAILED_READ_SEQUENCE_MATCH_X: return "failed_read_sequence_match_X";
        case VerifyResultType::FAILED_ASSEMBLY_SEAM_VALIDATION: return "failed_assembly_seam_validation";
        case VerifyResultType::FAILED_FINAL_SEQUENCE_MATCH: return "failed_final_sequence_match";
        case VerifyResultType::FAILED_MISSING_SEQUENCES: return "failed_missing_sequences";
        default: return "unknown";
    }
}

// Convert SegmentType to string
string segment_type_to_string(SegmentType type)
{
    switch (type) {
        case SegmentType::A: return "A";
        case SegmentType::SEAM_LEFT: return "seam_left";
        case SegmentType::X: return "X";
        case SegmentType::SEAM_RIGHT: return "seam_right";
        case SegmentType::B: return "B";
        default: return "unknown";
    }
}

// Helper function to format sequence info for debugging
string format_sequence_info(const string& seq, const string& label)
{
    if (seq.empty()) {
        return label + ": empty";
    }
    
    string first10 = seq.length() >= 10 ? seq.substr(0, 10) : seq;
    string last10 = seq.length() >= 10 ? seq.substr(seq.length() - 10) : seq;
    
    return label + ": len=" + to_string(seq.length()) + " first10=" + first10 + " last10=" + last10;
}

// Helper function to add a segment to the context with position-wise annotations
void add_segment(VerifyContext& ctx, SegmentType type, const string& sequence)
{
    if (!sequence.empty()) {
        uint32_t length_before = ctx.final_sequence.length();
        ctx.final_sequence += sequence;
        // add segment type annotation for each position in the added sequence
        ctx.segment_annotations.insert(ctx.segment_annotations.end(), sequence.length(), type);
        uint32_t length_after = ctx.final_sequence.length();
    }
}

// Helper function to remove a segment from the end of the context
void remove_segment(VerifyContext& ctx, uint32_t length)
{
    massert(length > 0 && length <= ctx.final_sequence.length(), 
            "remove_segment: invalid length=%u for final_sequence of length=%zu", 
            length, ctx.final_sequence.length());
    uint32_t length_before = ctx.final_sequence.length();
    ctx.final_sequence.resize(ctx.final_sequence.size() - length);
    ctx.segment_annotations.resize(ctx.segment_annotations.size() - length);
    uint32_t length_after = ctx.final_sequence.length();
}

// Helper function for detailed mismatch analysis
string analyze_sequence_mismatch(const string& constructed, const string& expected, const vector<SegmentType>& segment_annotations)
{
    size_t min_len = min(constructed.length(), expected.length());
    size_t mismatch_pos = min_len; // default to end if no mismatch found
    
    // find first mismatch position
    for (size_t i = 0; i < min_len; ++i) {
        if (constructed[i] != expected[i]) {
            mismatch_pos = i;
            break;
        }
    }
    
    // handle length mismatch case
    bool is_length_mismatch = (mismatch_pos == min_len && constructed.length() != expected.length());
    
    // find which segment contains the mismatch using position-wise annotations
    string segment_info = "unknown_segment";
    if (mismatch_pos < segment_annotations.size()) {
        SegmentType segment_type = segment_annotations[mismatch_pos];
        segment_info = segment_type_to_string(segment_type);
    } else if (is_length_mismatch && !segment_annotations.empty()) {
        // for length mismatches, report the last segment type
        SegmentType last_segment_type = segment_annotations.back();
        segment_info = "end_of_" + segment_type_to_string(last_segment_type);
    }
    
    // get context around mismatch (±30 nucleotides)
    size_t context_start = (mismatch_pos >= 30) ? mismatch_pos - 30 : 0;
    size_t context_end_constructed = min(mismatch_pos + 30, constructed.length());
    size_t context_end_expected = min(mismatch_pos + 30, expected.length());
    
    string constructed_context = constructed.substr(context_start, context_end_constructed - context_start);
    string expected_context = expected.substr(context_start, context_end_expected - context_start);
    
    // add length information for length mismatches
    string mismatch_description;
    if (is_length_mismatch) {
        mismatch_description = "length mismatch at position " + to_string(mismatch_pos) + 
                             " (constructed=" + to_string(constructed.length()) + 
                             ", expected=" + to_string(expected.length()) + ") in segment " + segment_info;
    } else {
        mismatch_description = "mismatch at position " + to_string(mismatch_pos) + " in segment " + segment_info;
    }
    
    return mismatch_description + 
           "\nconstructed: " + constructed_context + 
           "\nexpected:    " + expected_context;
}

RearrangeVerify::RearrangeVerify(const AssemblySequences* assembly_seqs, VerifyErrorMode mode)
    : assembly_sequences(assembly_seqs), error_mode(mode) {}


VerifyResult RearrangeVerify::verify_event(const ReadEvent& event, const string& lib_id, const AlignmentStore& store,
                                          const unordered_map<string, string>& sequences) const
{
    // basic validation - requires assembly sequences
    if (!can_verify()) {
        cerr << "error: verification requires assembly sequences to be loaded" << endl;
        return VerifyResult::failure(VerifyResultType::FAILED_MISSING_SEQUENCES, 
                                   "assembly sequences not loaded");
    }
    
    // check if required sequences exist
    if (!assembly_sequences->has_contig(event.contig_id)) {
        return VerifyResult::failure(VerifyResultType::FAILED_MISSING_SEQUENCES, 
                                   "missing contig sequence: " + event.contig_id);
    }
    
    

    // create verification context
    VerifyContext ctx;
    
    // step 1: create mutated sequences for A, B, and X
    VerifyResult result = create_mutated_sequences(event, store, ctx);
    if (result.is_failure()) {
        return result;
    }
    
    // step 2: verify read sequences match mutated assembly
    result = validate_read_sequences(event, lib_id, ctx, sequences);
    if (result.is_failure()) {
        return result;
    }
    
    // step 3: validate assembly seams (before clipping)
    result = validate_assembly_seams(event);
    if (result.is_failure()) {
        return result;
    }
    
    // steps 4-5: construct and validate final sequence
    result = construct_and_validate_final_sequence(event, ctx);
    if (result.is_failure()) {
        return result;
    }
    
    return VerifyResult::success();
}

VerifyResult RearrangeVerify::create_mutated_sequences(const ReadEvent& event, const AlignmentStore& store, VerifyContext& ctx) const
{
    // step 1: take assembly sequence of A, B, X and mutate based on mutations
    
    // get contig sequence
    string contig_sequence = assembly_sequences->get_sequence(event.contig_id);
    if (contig_sequence.empty()) {
        return VerifyResult::failure(VerifyResultType::FAILED_MISSING_SEQUENCES,
                                   "empty contig sequence for: " + event.contig_id);
    }
    
    // create A_tag by extracting and mutating alignment A
    const Alignment& A = event.alignment_A;
    string A_fragment = contig_sequence.substr(A.contig_start, A.contig_end - A.contig_start);
    ctx.A_tag = apply_mutations(A_fragment, A.mutations, store, A, event.read_id, event.contig_id);

    // create B_tag by extracting and mutating alignment B  
    const Alignment& B = event.alignment_B;
    string B_fragment = contig_sequence.substr(B.contig_start, B.contig_end - B.contig_start);
    ctx.B_tag = apply_mutations(B_fragment, B.mutations, store, B, event.read_id, event.contig_id);
 
    // create X_tag if present (for insertions and inversions)
    if (event.has_alignment_X) {
        const Alignment& X = event.alignment_X;
        string X_contig_sequence;
        
        // X might be on a different contig
        if (X.contig_index == A.contig_index) {
            X_contig_sequence = contig_sequence;
        } else {
            string X_contig_id = store.get_contig_id(X.contig_index);
            X_contig_sequence = assembly_sequences->get_sequence(X_contig_id);
            if (X_contig_sequence.empty()) {
                return VerifyResult::failure(VerifyResultType::FAILED_MISSING_SEQUENCES,
                                           "empty X contig sequence for: " + X_contig_id);
            }
        }
        
        string X_fragment = X_contig_sequence.substr(X.contig_start, X.contig_end - X.contig_start);
        ctx.X_tag = apply_mutations(X_fragment, X.mutations, store, X, event.read_id, 
                                   store.get_contig_id(X.contig_index));
    }
    
    return VerifyResult::success();
}

VerifyResult RearrangeVerify::validate_read_sequences(const ReadEvent& event, const string& /* lib_id */, VerifyContext& ctx,
                                                    const unordered_map<string, string>& sequences) const
{
    // step 2: verify read sequences of A_tag, B_tag and X_tag match the mutated assembly
    
    // get read sequence from passed sequences map
    auto seq_it = sequences.find(event.read_id);
    if (seq_it == sequences.end()) {
        return VerifyResult::failure(VerifyResultType::FAILED_MISSING_SEQUENCES,
                                   "read sequence not found: " + event.read_id);
    }
    
    ctx.read_sequence = seq_it->second;
    
    // validate A alignment
    const Alignment& A = event.alignment_A;
    string read_A = ctx.read_sequence.substr(A.read_start, A.read_end - A.read_start);
    if (A.is_reverse) {
        read_A = reverse_complement(read_A);
    }
    
    if (read_A != ctx.A_tag) {
        return VerifyResult::failure(VerifyResultType::FAILED_READ_SEQUENCE_MATCH_A,
                                   "A alignment sequence mismatch: " + format_sequence_info(read_A, "read_A") + 
                                   " vs " + format_sequence_info(ctx.A_tag, "A_tag"));
    }
    
    // validate B alignment  
    const Alignment& B = event.alignment_B;
    string read_B = ctx.read_sequence.substr(B.read_start, B.read_end - B.read_start);
    if (B.is_reverse) {
        read_B = reverse_complement(read_B);
    }
    
    if (read_B != ctx.B_tag) {
        return VerifyResult::failure(VerifyResultType::FAILED_READ_SEQUENCE_MATCH_B,
                                   "B alignment sequence mismatch: " + format_sequence_info(read_B, "read_B") + 
                                   " vs " + format_sequence_info(ctx.B_tag, "B_tag"));
    }
    
    // validate X alignment if present
    if (event.has_alignment_X) {
        const Alignment& X = event.alignment_X;
        string read_X = ctx.read_sequence.substr(X.read_start, X.read_end - X.read_start);
        if (X.is_reverse != (event.read_strand == "-")) {
            read_X = reverse_complement(read_X);
        }
        
        if (read_X != ctx.X_tag) {
            return VerifyResult::failure(VerifyResultType::FAILED_READ_SEQUENCE_MATCH_X,
                                       "X alignment sequence mismatch: " + format_sequence_info(read_X, "read_X") + 
                                       " vs " + format_sequence_info(ctx.X_tag, "X_tag"));
        }
    }
    
    return VerifyResult::success();
}

VerifyResult RearrangeVerify::validate_assembly_seams(const ReadEvent& event) const
{
    // step 3: validate assembly seams are in the assembly where we expect them
    
    if (event.assembly_seams.empty()) {
        return VerifyResult::success(); // no assembly seams to validate
    }
    
    // determine expected assembly seam count by event type
    int expected_assembly_seam_count;
    if (event.type == RearrangementType::LARGE_DELETE) {
        expected_assembly_seam_count = 0; // should be empty, but we already checked above
        return VerifyResult::success();
    } else if (event.type == RearrangementType::LARGE_INSERT) {
        expected_assembly_seam_count = 1;
    } else { // LARGE_INVERT
        expected_assembly_seam_count = 2;
    }
    
    // parse assembly seams
    vector<string> seams = parse_seams(event.assembly_seams, expected_assembly_seam_count);
    
    // get contig sequence
    string contig_sequence = assembly_sequences->get_sequence(event.contig_id);
    if (contig_sequence.empty()) {
        return VerifyResult::failure(VerifyResultType::FAILED_MISSING_SEQUENCES,
                                   "empty contig sequence for assembly seam validation: " + event.contig_id);
    }
    
    // for insertions: assembly seam should be between clipped A and B coordinates
    // for inversions: assembly seams should be A:X and X:B intervals
    // for deletions: no assembly seams expected
    
    if (event.type == RearrangementType::LARGE_INSERT) {
        // assembly seam should be the sequence between out_clip and in_clip
        if (seams.size() != 1) {
            return VerifyResult::failure(VerifyResultType::FAILED_ASSEMBLY_SEAM_VALIDATION,
                                       "insertion should have exactly 1 assembly seam, found: " + to_string(seams.size()));
        }
        
        // assert that in_clip >= out_clip for valid substring
        massert(event.in_clip >= event.out_clip, "in_clip (%u) must be >= out_clip (%u)", event.in_clip, event.out_clip);
        
        string expected_seam = contig_sequence.substr(event.out_clip, event.in_clip - event.out_clip);
        string actual_seam = extract_seam_sequence(seams[0]);
        
        if (actual_seam != expected_seam) {
            return VerifyResult::failure(VerifyResultType::FAILED_ASSEMBLY_SEAM_VALIDATION,
                                       "insertion assembly seam mismatch: " + format_sequence_info(actual_seam, "actual") + 
                                       " vs " + format_sequence_info(expected_seam, "expected"));
        }
    }
    else if (event.type == RearrangementType::LARGE_INVERT) {
        // should have two seams: A:X and X:B
        if (seams.size() != 2) {
            return VerifyResult::failure(VerifyResultType::FAILED_ASSEMBLY_SEAM_VALIDATION,
                                       "inversion should have exactly 2 assembly seams, found: " + to_string(seams.size()));
        }
        
        // validate the seam sequences exist in the expected locations
        // this is a simplified validation - more complex logic would verify exact positions
        for (const string& seam : seams) {
            string seam_seq = extract_seam_sequence(seam);
            if (!seam_seq.empty() && contig_sequence.find(seam_seq) == string::npos) {
                return VerifyResult::failure(VerifyResultType::FAILED_ASSEMBLY_SEAM_VALIDATION,
                                           "inversion seam sequence not found in contig: " + format_sequence_info(seam_seq, "seam"));
            }
        }
    }
    
    return VerifyResult::success();
}

VerifyResult RearrangeVerify::construct_and_validate_final_sequence(const ReadEvent& event, VerifyContext& ctx) const
{
    const Alignment& A = event.alignment_A;
    const Alignment& B = event.alignment_B;
    
    // determine expected seam count by event type
    int expected_seam_count;
    if (event.type == RearrangementType::LARGE_DELETE) {
        expected_seam_count = 1;
    } else {
        expected_seam_count = 2; // insert/invert events
    }
    
    // parse read seams into structured format with validation
    vector<Seam> read_seams = parse_seams_structured(event.read_seams, expected_seam_count);
    
    // construct sequence based on event type using position-wise annotations
    ctx.final_sequence.clear();
    ctx.segment_annotations.clear();
    
    if (event.type == RearrangementType::LARGE_DELETE) {
        // deletion: A_tag + read_seam + B_tag
        add_segment(ctx, SegmentType::A, ctx.A_tag);
        
        // handle seam (gap or overlap)
        const Seam& seam = read_seams[0];
        if (seam.is_gap && !seam.sequence.empty()) {
            add_segment(ctx, SegmentType::SEAM_LEFT, seam.sequence);
        } else if (!seam.sequence.empty()) {
            // overlap: verify exact match and remove from end of sequence
            massert(seam.sequence.length() <= ctx.final_sequence.length(),
                    "overlap seam length (%zu) exceeds current sequence length (%zu)",
                    seam.sequence.length(), ctx.final_sequence.length());
            string actual_suffix = ctx.final_sequence.substr(ctx.final_sequence.size() - seam.sequence.length());
            if (actual_suffix == seam.sequence) {
                remove_segment(ctx, seam.sequence.length());
            } else {
                return VerifyResult::failure(VerifyResultType::FAILED_FINAL_SEQUENCE_MATCH,
                                             "overlap mismatch in delete seam: expected=" + seam.sequence + " actual=" + actual_suffix);
            }
        }
        
        add_segment(ctx, SegmentType::B, ctx.B_tag);
    }
    else if (event.type == RearrangementType::LARGE_INSERT) {
        // insertion: A_tag + read_seam_left + X_tag + read_seam_right + B_tag
        add_segment(ctx, SegmentType::A, ctx.A_tag);
        
        // handle left seam (gap or overlap)
        const Seam& left_seam = read_seams[0];
        if (left_seam.is_gap && !left_seam.sequence.empty()) {
            add_segment(ctx, SegmentType::SEAM_LEFT, left_seam.sequence);
        } else if (!left_seam.sequence.empty()) {
            // overlap: verify exact match and remove from end of sequence
            if (left_seam.sequence.length() <= ctx.final_sequence.length()) {
                string actual_suffix = ctx.final_sequence.substr(ctx.final_sequence.size() - left_seam.sequence.length());
                if (actual_suffix == left_seam.sequence) {
                    remove_segment(ctx, left_seam.sequence.length());
                } else {
                    cout << "overlap mismatch in insert left seam: expected=" << left_seam.sequence << " actual=" << actual_suffix << endl;
                    return VerifyResult::failure(VerifyResultType::FAILED_FINAL_SEQUENCE_MATCH, "overlap seam mismatch");
                }
            }
        }
        
        // apply reverse complement to X_tag if X alignment is on negative strand
        string X_final = ctx.X_tag;
        if (event.has_alignment_X && event.alignment_X.is_reverse) {
            X_final = reverse_complement(X_final);
        }        
        add_segment(ctx, SegmentType::X, X_final);
        
        // handle right seam (gap or overlap)
        const Seam& right_seam = read_seams[1];
        if (right_seam.is_gap && !right_seam.sequence.empty()) {
            add_segment(ctx, SegmentType::SEAM_RIGHT, right_seam.sequence);
        } else if (!right_seam.sequence.empty()) {
            // overlap: verify exact match and remove from end of sequence
            if (right_seam.sequence.length() <= ctx.final_sequence.length()) {
                string actual_suffix = ctx.final_sequence.substr(ctx.final_sequence.size() - right_seam.sequence.length());
                if (actual_suffix == right_seam.sequence) {
                    remove_segment(ctx, right_seam.sequence.length());
                } else {
                    cout << "overlap mismatch in insert right seam: expected=" << right_seam.sequence << " actual=" << actual_suffix << endl;
                    return VerifyResult::failure(VerifyResultType::FAILED_FINAL_SEQUENCE_MATCH, "overlap seam mismatch");
                }
            }
        }
        
        add_segment(ctx, SegmentType::B, ctx.B_tag);
    }
    else if (event.type == RearrangementType::LARGE_INVERT) {
        // inversion: A_tag + read_seam_left + inverted_X_tag + read_seam_right + B_tag
        add_segment(ctx, SegmentType::A, ctx.A_tag);

        // handle left seam (gap or overlap)
        const Seam& left_seam = read_seams[0];
        if (left_seam.is_gap && !left_seam.sequence.empty()) {
            add_segment(ctx, SegmentType::SEAM_LEFT, left_seam.sequence);
        } else if (!left_seam.sequence.empty()) {
            // overlap: verify exact match and remove from end of sequence
            if (left_seam.sequence.length() <= ctx.final_sequence.length()) {
                string actual_suffix = ctx.final_sequence.substr(ctx.final_sequence.size() - left_seam.sequence.length());
                if (actual_suffix == left_seam.sequence) {
                    remove_segment(ctx, left_seam.sequence.length());
                } else {
                    cout << "overlap mismatch in invert left seam: expected=" << left_seam.sequence << " actual=" << actual_suffix << endl;
                    return VerifyResult::failure(VerifyResultType::FAILED_FINAL_SEQUENCE_MATCH, "overlap seam mismatch");
                }
            }
        }
        
        // apply reverse complement to X_tag if X alignment is on negative strand
        string X_final = ctx.X_tag;
        if (event.has_alignment_X && event.alignment_X.is_reverse) {
            X_final = reverse_complement(X_final);
        }
        add_segment(ctx, SegmentType::X, X_final);
        
        // handle right seam (gap or overlap)
        const Seam& right_seam = read_seams[1];
        if (right_seam.is_gap && !right_seam.sequence.empty()) {
            add_segment(ctx, SegmentType::SEAM_RIGHT, right_seam.sequence);
        } else if (!right_seam.sequence.empty()) {
            // overlap: verify exact match and remove from end of sequence
            if (right_seam.sequence.length() <= ctx.final_sequence.length()) {
                string actual_suffix = ctx.final_sequence.substr(ctx.final_sequence.size() - right_seam.sequence.length());
                if (actual_suffix == right_seam.sequence) {
                    remove_segment(ctx, right_seam.sequence.length());
                } else {
                    cout << "overlap mismatch in invert right seam: expected=" << right_seam.sequence << " actual=" << actual_suffix << endl;
                    return VerifyResult::failure(VerifyResultType::FAILED_FINAL_SEQUENCE_MATCH, "overlap seam mismatch");
                }
            }
        }
        
        add_segment(ctx, SegmentType::B, ctx.B_tag);
    }
    
    // step 5: validate final sequence against read sequence
    // extract the relevant portion of the read sequence
    uint32_t read_start, read_end;
    bool is_read_reversed = (event.read_strand == "-");
    
    if (is_read_reversed) {
        // for reversed reads: sequence goes from start of B to end of A
        read_start = B.read_start;
        read_end = A.read_end;
    } else {
        // for forward reads: sequence goes from start of A to end of B
        read_start = A.read_start;
        read_end = B.read_end;
    }
    
    if (read_end > ctx.read_sequence.length() || read_start >= read_end) {
        return VerifyResult::failure(VerifyResultType::FAILED_FINAL_SEQUENCE_MATCH,
                                   "invalid read coordinates: start=" + to_string(read_start) + 
                                   " end=" + to_string(read_end) + " read_len=" + to_string(ctx.read_sequence.length()) +
                                   " is_read_reversed=" + (is_read_reversed ? "true" : "false"));
    }
    
    string expected_read_segment = ctx.read_sequence.substr(read_start, read_end - read_start);
    
    // apply reverse complement if read is on reverse strand
    if (is_read_reversed) {
        expected_read_segment = reverse_complement(expected_read_segment);
    }
    
    if (ctx.final_sequence != expected_read_segment) {
        string detailed_analysis = analyze_sequence_mismatch(ctx.final_sequence, expected_read_segment, ctx.segment_annotations);
        string segment_report = create_segment_length_report(event, ctx);
        return VerifyResult::failure(VerifyResultType::FAILED_FINAL_SEQUENCE_MATCH,
                                   "final sequence mismatch: " + format_sequence_info(ctx.final_sequence, "constructed") + 
                                   " vs " + format_sequence_info(expected_read_segment, "expected=" + to_string(expected_read_segment.length())) + 
                                   " | " + segment_report + " | " + detailed_analysis);
    }
    
    return VerifyResult::success();
}

void RearrangeVerify::verify_all_events(const string& lib_id, vector<ReadEvent>& events, const AlignmentStore& store,
                                        const unordered_map<string, string>& sequences)
{
    cout << "verifying " << events.size() << " events..." << endl;
    
    size_t verified_events = 0;
    size_t failed_events = 0;
    size_t detailed_reports = 0;
    const size_t max_detailed_reports = 10;
    
    // iterate through events and verify each one
    size_t event_index = 0;
    auto it = events.begin();
    while (it != events.end()) {
        // RearrangeReadEvent::print_event(*it, store);

        VerifyResult result = verify_event(*it, lib_id, store, sequences);

        if (result.is_success()) {
            verified_events++;
            ++it; // keep the event
        } else {
            failed_events++;

            // print detailed report for first 3 failures
            if (detailed_reports < max_detailed_reports) {
                cout << "==============================================" << endl;
                cout << "VERIFY FAILED [" << (detailed_reports + 1) << "/3]: index=" << event_index
                     << " read_id=" << it->read_id 
                     << " event_type=" << it->type 
                     << " reason=" << verify_result_type_to_string(result.type) << endl;
                cout << "  message: " << result.message << endl;
                RearrangeReadEvent::print_event(*it, store);
                detailed_reports++;
            }

            // handle error modes
            if (error_mode == VerifyErrorMode::EXIT_ON_ERROR) {
                cout << "verification failed, exiting on first error" << endl;
                exit(1);
            } else if (error_mode == VerifyErrorMode::EXIT_IF_ERROR) {
                // continue processing but will exit at the end
            } else { // WARNING_ONLY
                // just continue
            }
        }
        ++event_index;
    }
    
    cout << "verification complete: " << verified_events << " passed, " << failed_events << " failed" << endl;
    
    // exit if error mode requires it and we had failures
    if (failed_events > 0 && error_mode == VerifyErrorMode::EXIT_IF_ERROR) {
        cout << "verification had " << failed_events << " failures, exiting" << endl;
        exit(1);
    }
}

vector<string> RearrangeVerify::parse_seams(const string& seam_string, int expected_seam_count) const
{
    vector<string> seams;
    
    // count colons to validate expected seam count
    int colon_count = 0;
    for (char c : seam_string) {
        if (c == ':') colon_count++;
    }
    
    // validate colon count matches expected seam count
    int expected_colons = expected_seam_count - 1;
    massert(colon_count == expected_colons, 
           "seam string has %d colons, expected %d for %d seams: '%s'", 
           colon_count, expected_colons, expected_seam_count, seam_string.c_str());
    
    if (seam_string.empty()) {
        // empty string should only be valid for 1 expected seam
        massert(expected_seam_count == 1, "empty seam string only valid for 1 expected seam, got %d", expected_seam_count);
        seams.push_back("");
        return seams;
    }
    
    // split by ':'
    size_t start = 0;
    size_t pos = 0;
    while ((pos = seam_string.find(':', start)) != string::npos) {
        seams.push_back(seam_string.substr(start, pos - start));
        start = pos + 1;
    }
    seams.push_back(seam_string.substr(start)); // last part
    
    // ensure we have exactly the expected number of seams
    massert(seams.size() == static_cast<size_t>(expected_seam_count), 
           "parsed %zu seams, expected %d from string: '%s'", 
           seams.size(), expected_seam_count, seam_string.c_str());
    
    return seams;
}

vector<Seam> RearrangeVerify::parse_seams_structured(const string& seam_string, int expected_seam_count) const
{
    // first parse as strings with validation
    vector<string> seam_strings = parse_seams(seam_string, expected_seam_count);
    
    // convert to structured seams
    vector<Seam> seams;
    for (const string& seam_str : seam_strings) {
        bool is_gap = is_gap_seam(seam_str);
        string sequence = extract_seam_sequence(seam_str);
        seams.emplace_back(is_gap, sequence);
    }
    
    return seams;
}

bool RearrangeVerify::is_gap_seam(const string& seam) const
{
    if (seam.empty()) return false;
    if (seam[0] == '+') return true;   // gap
    if (seam[0] == '-') return false;  // overlap
    return false;  // no prefix means empty seam
}

string RearrangeVerify::extract_seam_sequence(const string& seam) const
{
    if (seam.empty()) {
        return seam; // empty seam
    }
    if (seam[0] == '+' || seam[0] == '-') {
        return seam.substr(1); // remove prefix for gaps and overlaps
    }
    return seam; // no prefix means empty seam
}

int32_t RearrangeVerify::calculate_seam_read_length(const Alignment& align1, const Alignment& align2, bool is_read_reversed) const
{
    int32_t result;
    if (is_read_reversed) {
        // for reversed reads, order is flipped (B, X, A), so swap the alignments
        result = static_cast<int32_t>(align1.read_start) - static_cast<int32_t>(align2.read_end);
    } else {
        // for forward reads, normal order (A, X, B)
        result = static_cast<int32_t>(align2.read_start) - static_cast<int32_t>(align1.read_end);
    }
    
    massert(result >= 0, "seam read length is negative (%d): align1(idx=%u %u-%u) align2(idx=%u %u-%u) is_read_reversed=%s",
           result, align1.read_index, align1.read_start, align1.read_end,
           align2.read_index, align2.read_start, align2.read_end,
           is_read_reversed ? "true" : "false");
    
    return result;
}

string RearrangeVerify::create_segment_length_report(const ReadEvent& event, const VerifyContext& ctx) const
{
    const Alignment& A = event.alignment_A;
    const Alignment& B = event.alignment_B;
    
    // determine if read is reversed
    bool is_read_reversed = (event.read_strand == "-");
    
    // determine expected seam count by event type
    int expected_seam_count;
    if (event.type == RearrangementType::LARGE_DELETE) {
        expected_seam_count = 1;
    } else {
        expected_seam_count = 2; // insert/invert events
    }
    
    // parse read seams to get gap/overlap info
    vector<Seam> read_seams = parse_seams_structured(event.read_seams, expected_seam_count);
    
    vector<string> segment_reports;
    
    // A segment: C=constructed length, R=read length
    uint32_t A_constructed = ctx.A_tag.length();
    uint32_t A_read = A.read_end - A.read_start;
    segment_reports.push_back("A=(C=" + to_string(A_constructed) + ",R=" + to_string(A_read) + ")");
    
    // seam segments based on event type
    if (event.type == RearrangementType::LARGE_DELETE) {
        if (!read_seams.empty()) {
            uint32_t seam_constructed = read_seams[0].sequence.length();
            // read length is the gap/overlap between A and B
            int32_t seam_read = calculate_seam_read_length(A, B, is_read_reversed);
            string seam_type = read_seams[0].is_gap ? "gap" : "overlap";
            segment_reports.push_back("SEAM=(" + seam_type + "," + to_string(seam_constructed) + "," + to_string(seam_read) + ")");
        }
    }
    else if (event.type == RearrangementType::LARGE_INSERT || event.type == RearrangementType::LARGE_INVERT) {
        const Alignment& X = event.alignment_X;
        
        // left seam
        if (read_seams.size() >= 1) {
            uint32_t seam_constructed = read_seams[0].sequence.length();
            int32_t seam_read = calculate_seam_read_length(A, X, is_read_reversed);
            string seam_type = read_seams[0].is_gap ? "gap" : "overlap";
            segment_reports.push_back("SEAM_LEFT=(" + seam_type + "," + to_string(seam_constructed) + "," + to_string(seam_read) + ")");
        }
        
        // X segment
        uint32_t X_constructed = ctx.X_tag.length();
        uint32_t X_read = X.read_end - X.read_start;
        segment_reports.push_back("X=(C=" + to_string(X_constructed) + ",R=" + to_string(X_read) + ")");
        
        // right seam
        if (read_seams.size() >= 2) {
            uint32_t seam_constructed = read_seams[1].sequence.length();
            int32_t seam_read = calculate_seam_read_length(X, B, is_read_reversed);
            string seam_type = read_seams[1].is_gap ? "gap" : "overlap";
            segment_reports.push_back("SEAM_RIGHT=(" + seam_type + "," + to_string(seam_constructed) + "," + to_string(seam_read) + ")");
        }
    }
    
    // B segment
    uint32_t B_constructed = ctx.B_tag.length();
    uint32_t B_read = B.read_end - B.read_start;
    segment_reports.push_back("B=(C=" + to_string(B_constructed) + ",R=" + to_string(B_read) + ")");
    
    // join all segments
    string result = "Segment lengths: ";
    for (size_t i = 0; i < segment_reports.size(); ++i) {
        if (i > 0) result += " ";
        result += segment_reports[i];
    }
    
    // add total lengths
    result += " | Total: constructed=" + to_string(ctx.final_sequence.length());
    
    return result;
}

