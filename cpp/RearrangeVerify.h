#ifndef REARRANGE_VERIFY_H
#define REARRANGE_VERIFY_H

#include <string>
#include <unordered_map>
#include "Sequences.h"

// Forward declarations
struct ReadEvent;
struct Alignment;
class AlignmentStore;

// Segment types for tracking constructed sequence components
enum class SegmentType {
    A,
    SEAM_LEFT,
    X,
    SEAM_RIGHT,
    B
};

// Segment information for mismatch reporting
struct ReadSegment {
    SegmentType type;
    uint32_t start_pos;  // start position in final sequence
    uint32_t length;     // length of this segment
    
    ReadSegment(SegmentType t, uint32_t start, uint32_t len) : type(t), start_pos(start), length(len) {}
};

// Convert SegmentType to string
std::string segment_type_to_string(SegmentType type);

// Seam structure to distinguish between gaps and overlaps
struct Seam {
    bool is_gap;        // true for gaps (+sequence), false for overlaps (sequence)
    std::string sequence; // the actual sequence (without + prefix)
    
    Seam(bool gap, const std::string& seq) : is_gap(gap), sequence(seq) {}
};

// Verification context to hold intermediate results
struct VerifyContext {
    std::string A_tag;
    std::string B_tag;
    std::string X_tag;
    std::string read_sequence;
    std::string final_sequence;
    std::vector<SegmentType> segment_annotations;  // position-wise segment annotations, same length as final_sequence
};

// Verification result types for detailed failure tracking
enum class VerifyResultType {
    SUCCESS,
    FAILED_ASSEMBLY_MUTATION_A,
    FAILED_ASSEMBLY_MUTATION_B, 
    FAILED_ASSEMBLY_MUTATION_X,
    FAILED_READ_SEQUENCE_MATCH_A,
    FAILED_READ_SEQUENCE_MATCH_B,
    FAILED_READ_SEQUENCE_MATCH_X,
    FAILED_ASSEMBLY_SEAM_VALIDATION,
    FAILED_FINAL_SEQUENCE_MATCH,
    FAILED_MISSING_SEQUENCES
};

// Verification result class with type and detailed message
class VerifyResult {
public:
    VerifyResultType type;
    std::string message;
    
    VerifyResult(VerifyResultType t, const std::string& msg = "") : type(t), message(msg) {}
    
    // convenience constructors
    static VerifyResult success() { return VerifyResult(VerifyResultType::SUCCESS, "verification successful"); }
    static VerifyResult failure(VerifyResultType t, const std::string& msg) { return VerifyResult(t, msg); }
    
    // check if successful
    bool is_success() const { return type == VerifyResultType::SUCCESS; }
    bool is_failure() const { return type != VerifyResultType::SUCCESS; }
};


// Verification error handling modes
enum class VerifyErrorMode {
    EXIT_ON_ERROR,    // stop with clear messages on first error
    EXIT_IF_ERROR,    // complete all events, report stats, exit if any errors
    WARNING_ONLY      // does not break if errors were found
};

// Conversion functions for VerifyErrorMode
VerifyErrorMode string_to_verify_error_mode(const std::string& str);
std::string verify_error_mode_to_string(VerifyErrorMode mode);

// Conversion function for VerifyResultType
std::string verify_result_type_to_string(VerifyResultType result_type);

class RearrangeVerify {
private:
    const AssemblySequences* assembly_sequences;
    VerifyErrorMode error_mode;

public:
    RearrangeVerify(const AssemblySequences* assembly_seqs = nullptr, 
                   VerifyErrorMode mode = VerifyErrorMode::EXIT_ON_ERROR);
    
    // verify a read event
    VerifyResult verify_event(const ReadEvent& event, const std::string& lib_id, const AlignmentStore& store,
                             const std::unordered_map<std::string, std::string>& sequences) const;
    
    // verify all events after seam resolution
    void verify_all_events(const std::string& lib_id, std::vector<ReadEvent>& events, const AlignmentStore& store,
                           const std::unordered_map<std::string, std::string>& sequences);
    
    // accessors
    VerifyErrorMode get_error_mode() const { return error_mode; }
    
    // check if verification is possible
    bool can_verify() const { return assembly_sequences && assembly_sequences->is_loaded(); }
    
private:
    // step-by-step verification functions
    VerifyResult create_mutated_sequences(const ReadEvent& event, const AlignmentStore& store, VerifyContext& ctx) const;
    VerifyResult validate_read_sequences(const ReadEvent& event, const std::string& lib_id, VerifyContext& ctx,
                                        const std::unordered_map<std::string, std::string>& sequences) const;
    VerifyResult validate_assembly_seams(const ReadEvent& event) const;
    VerifyResult construct_and_validate_final_sequence(const ReadEvent& event, VerifyContext& ctx) const;
    
    // seam parsing helpers
    std::vector<std::string> parse_seams(const std::string& seam_string, int expected_seam_count) const;
    std::vector<Seam> parse_seams_structured(const std::string& seam_string, int expected_seam_count) const;
    bool is_gap_seam(const std::string& seam) const;  // starts with '+' (gap) vs '-' (overlap)
    std::string extract_seam_sequence(const std::string& seam) const;
    
    // debug helpers
    std::string create_segment_length_report(const ReadEvent& event, const VerifyContext& ctx) const;
    int32_t calculate_seam_read_length(const Alignment& align1, const Alignment& align2, bool is_read_reversed) const;
};

#endif // REARRANGE_VERIFY_H
