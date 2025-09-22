#include "RearrangeErrors.h"
#include <iostream>

using namespace std;

RearrangeErrors::RearrangeErrors()
{
    // initialize rejection counts to zero
    rejection_counts[EventTestResult::REJECT_ALIGNMENT_OVERLAP] = 0;
    rejection_counts[EventTestResult::REJECT_INSERT_ELEMENT_REINSERTED] = 0;
    rejection_counts[EventTestResult::REJECT_INSERT_ELEMENT_OVERLAPS] = 0;
    rejection_counts[EventTestResult::REJECT_ELEMENT_TOO_SHORT] = 0;
}

void RearrangeErrors::record_rejection(const std::string& read_id, EventTestResult reason,
                                      const Alignment& A, const Alignment& B, 
                                      const AlignmentStore& store,
                                      const Alignment* X)
{
    RejectionRecord record(read_id, reason);
    
    // record alignment A details
    record.A_contig_start = A.contig_start;
    record.A_contig_end = A.contig_end;
    record.A_read_start = A.read_start;
    record.A_read_end = A.read_end;
    record.A_is_reverse = A.is_reverse;
    record.A_contig_id = store.get_contig_id(A.contig_index);
    
    // record alignment B details
    record.B_contig_start = B.contig_start;
    record.B_contig_end = B.contig_end;
    record.B_read_start = B.read_start;
    record.B_read_end = B.read_end;
    record.B_is_reverse = B.is_reverse;
    record.B_contig_id = store.get_contig_id(B.contig_index);
    
    // record alignment X details if present
    if (X) {
        record.has_X = true;
        record.X_contig_start = X->contig_start;
        record.X_contig_end = X->contig_end;
        record.X_read_start = X->read_start;
        record.X_read_end = X->read_end;
        record.X_is_reverse = X->is_reverse;
        record.X_contig_id = store.get_contig_id(X->contig_index);
    }
    
    rejection_records.push_back(record);
    rejection_counts[reason]++;
}

void RearrangeErrors::write_to_file(const std::string& lib_prefix)
{
    string filename = lib_prefix + "_rejection_details.tsv";
    cout << "writing rejection details to " << filename << endl;
    
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "error: could not open file " << filename << endl;
        return;
    }
    
    // write header
    file << "read_id\trejection_reason\t"
         << "A_contig_id\tA_contig_start\tA_contig_end\tA_read_start\tA_read_end\tA_is_reverse\t"
         << "B_contig_id\tB_contig_start\tB_contig_end\tB_read_start\tB_read_end\tB_is_reverse\t"
         << "has_X\tX_contig_id\tX_contig_start\tX_contig_end\tX_read_start\tX_read_end\tX_is_reverse\n";
    
    // write data rows
    for (const RejectionRecord& record : rejection_records) {
        file << record.read_id << "\t"
             << record.rejection_reason << "\t"
             << record.A_contig_id << "\t"
             << record.A_contig_start << "\t"
             << record.A_contig_end << "\t"
             << record.A_read_start << "\t"
             << record.A_read_end << "\t"
             << (record.A_is_reverse ? "true" : "false") << "\t"
             << record.B_contig_id << "\t"
             << record.B_contig_start << "\t"
             << record.B_contig_end << "\t"
             << record.B_read_start << "\t"
             << record.B_read_end << "\t"
             << (record.B_is_reverse ? "true" : "false") << "\t"
             << (record.has_X ? "true" : "false") << "\t";
        
        if (record.has_X) {
            file << record.X_contig_id << "\t"
                 << record.X_contig_start << "\t"
                 << record.X_contig_end << "\t"
                 << record.X_read_start << "\t"
                 << record.X_read_end << "\t"
                 << (record.X_is_reverse ? "true" : "false");
        } else {
            file << "NA\t0\t0\t0\t0\tfalse";
        }
        file << "\n";
    }
    
    file.close();
    cout << "wrote " << rejection_records.size() << " rejection records to " << filename << endl;
}

void RearrangeErrors::print_summary() const
{
    cout << "rejection summary:" << endl;
    for (const auto& pair : rejection_counts) {
        cout << "  " << pair.first << ": " << pair.second << endl;
    }
    cout << "  total rejections: " << get_total_rejections() << endl;
}

size_t RearrangeErrors::get_rejection_count(EventTestResult reason) const
{
    auto it = rejection_counts.find(reason);
    return (it != rejection_counts.end()) ? it->second : 0;
}

size_t RearrangeErrors::get_total_rejections() const
{
    size_t total = 0;
    for (const auto& pair : rejection_counts) {
        total += pair.second;
    }
    return total;
}

const std::vector<RejectionRecord>& RearrangeErrors::get_rejection_records() const
{
    return rejection_records;
}

const std::map<EventTestResult, size_t>& RearrangeErrors::get_rejection_counts() const
{
    return rejection_counts;
}
