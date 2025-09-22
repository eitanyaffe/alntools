#ifndef REARRANGE_ERRORS_H
#define REARRANGE_ERRORS_H

#include "aln_types.h"
#include "alignment_store.h"
#include <string>
#include <vector>
#include <fstream>
#include <map>

enum class EventTestResult {
    FOUND_INSERTION,
    FOUND_INVERSION, 
    FOUND_DELETION,
    REJECT_ALIGNMENT_OVERLAP,
    REJECT_INSERT_ELEMENT_REINSERTED,
    REJECT_INSERT_ELEMENT_OVERLAPS,
    REJECT_ELEMENT_TOO_SHORT,
    REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE,
    REJECT_READ_GAP_TOO_LARGE
};

// add operator<< for EventTestResult
inline std::ostream& operator<<(std::ostream& os, const EventTestResult& result)
{
    switch (result) {
    case EventTestResult::FOUND_INSERTION:
        os << "found_insertion";
        break;
    case EventTestResult::FOUND_INVERSION:
        os << "found_inversion";
        break;
    case EventTestResult::FOUND_DELETION:
        os << "found_deletion";
        break;
    case EventTestResult::REJECT_ALIGNMENT_OVERLAP:
        os << "reject_alignment_overlap";
        break;
    case EventTestResult::REJECT_INSERT_ELEMENT_REINSERTED:
        os << "reject_insert_element_reinserted";
        break;
    case EventTestResult::REJECT_INSERT_ELEMENT_OVERLAPS:
        os << "reject_insert_element_overlaps";
        break;
    case EventTestResult::REJECT_ELEMENT_TOO_SHORT:
        os << "reject_element_too_short";
        break;
    case EventTestResult::REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE:
        os << "reject_contained_element_margin_too_large";
        break;
    case EventTestResult::REJECT_READ_GAP_TOO_LARGE:
        os << "reject_read_gap_too_large";
        break;
    default:
        os << "unknown";
        break;
    }
    return os;
}

struct RejectionRecord {
    std::string read_id;
    EventTestResult rejection_reason;
    
    // alignment A details
    uint32_t A_contig_start;
    uint32_t A_contig_end;
    uint32_t A_read_start;
    uint32_t A_read_end;
    bool A_is_reverse;
    std::string A_contig_id;
    
    // alignment B details
    uint32_t B_contig_start;
    uint32_t B_contig_end;
    uint32_t B_read_start;
    uint32_t B_read_end;
    bool B_is_reverse;
    std::string B_contig_id;
    
    // alignment X details (optional)
    bool has_X;
    uint32_t X_contig_start;
    uint32_t X_contig_end;
    uint32_t X_read_start;
    uint32_t X_read_end;
    bool X_is_reverse;
    std::string X_contig_id;
    
    RejectionRecord(const std::string& read_id = "", EventTestResult reason = EventTestResult::REJECT_ALIGNMENT_OVERLAP)
        : read_id(read_id), rejection_reason(reason), has_X(false) {}
};

class RearrangeErrors {
private:
    std::vector<RejectionRecord> rejection_records;
    std::map<EventTestResult, size_t> rejection_counts;

public:
    RearrangeErrors();
    
    void record_rejection(const std::string& read_id, EventTestResult reason,
                         const Alignment& A, const Alignment& B, 
                         const AlignmentStore& store,
                         const Alignment* X = nullptr);
    
    void write_to_file(const std::string& lib_prefix);
    
    void print_summary() const;
    
    size_t get_rejection_count(EventTestResult reason) const;
    size_t get_total_rejections() const;
    
    // accessors for R interface
    const std::vector<RejectionRecord>& get_rejection_records() const;
    const std::map<EventTestResult, size_t>& get_rejection_counts() const;
};

#endif // REARRANGE_ERRORS_H
