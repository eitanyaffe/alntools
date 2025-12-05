#include "Segmentation.h"
#include <algorithm>
#include <iostream>
#include <set>

using namespace std;

Segmentation::Segmentation(const AlignmentStore& store,
                         int min_anchor_length,
                         int min_dangle_length,
                         double max_anchor_mutations_percent,
                         int min_alignment_distance)
    : store(store)
    , min_anchor_length(min_anchor_length)
    , min_dangle_length(min_dangle_length)
    , max_anchor_mutations_percent(max_anchor_mutations_percent)
    , min_alignment_distance(min_alignment_distance)
{
}

bool Segmentation::is_valid_anchor(const Alignment& aln) const
{
    uint32_t alignment_length = aln.read_end - aln.read_start;
    if (static_cast<int>(alignment_length) < min_anchor_length) {
        return false;
    }
    
    double mutation_percent = (static_cast<double>(aln.get_mutation_count()) / 
        static_cast<double>(alignment_length)) * 100.0;
    if (mutation_percent > max_anchor_mutations_percent) {
        return false;
    }
    
    return true;
}

string Segmentation::get_breakpoint_type(bool is_reverse, bool is_end) const
{
    if (is_reverse) {
        return is_end ? "dangle_left" : "dangle_right";
    } else {
        return is_end ? "dangle_right" : "dangle_left";
    }
}

uint32_t Segmentation::get_breakpoint_coord(const Alignment& aln, bool is_end) const
{
    return is_end ? aln.contig_end : aln.contig_start;
}

bool Segmentation::is_too_close_to_read_boundary(const Alignment& anchor, uint32_t read_length,
						 bool check_next) const
{
    int threshold = min_alignment_distance + min_dangle_length;
    
    if (check_next) {
        uint32_t distance_to_end = read_length - anchor.read_end;
        return static_cast<int>(distance_to_end) < threshold;
    } else {
        uint32_t distance_to_start = anchor.read_start;
        return static_cast<int>(distance_to_start) < threshold;
    }
}

bool Segmentation::has_nearby_same_contig_alignment(const Alignment& anchor,
                                                     size_t anchor_idx,
                                                     bool check_next,
                                                     const vector<size_t>& alignment_indices,
                                                     const vector<Alignment>& alignments) const
{
    int threshold_read_distance = min_alignment_distance;
    uint32_t anchor_boundary = check_next ? anchor.read_end : anchor.read_start;
    
    int start_idx = check_next ? (anchor_idx + 1) : (anchor_idx - 1);
    
    if (start_idx < 0 || start_idx >= static_cast<int>(alignment_indices.size())) {
        return false;
    }
    
    int step = check_next ? 1 : -1;
    
    for (int idx = start_idx; idx >= 0 && idx < static_cast<int>(alignment_indices.size()); idx += step) {
        const Alignment& candidate = alignments[alignment_indices[idx]];
        
        int read_distance;
        if (check_next) {
            read_distance = static_cast<int>(candidate.read_start) - static_cast<int>(anchor_boundary);
        } else {
            read_distance = static_cast<int>(anchor_boundary) - static_cast<int>(candidate.read_end);
        }
        
        if (read_distance < 0) {
            continue;
        }
        
        if (read_distance > threshold_read_distance) {
            break;
        }
        
        if (anchor.contig_index != candidate.contig_index || anchor.is_reverse != candidate.is_reverse) {
            continue;
        }
        
        uint32_t candidate_length = candidate.read_end - candidate.read_start;
        if (static_cast<int>(candidate_length) < min_dangle_length) {
            continue;
        }
        
        uint32_t contig_distance = 0;
        if (check_next) {
            if (anchor.contig_end <= candidate.contig_start) {
                contig_distance = candidate.contig_start - anchor.contig_end;
            }
        } else {
            if (candidate.contig_end <= anchor.contig_start) {
                contig_distance = anchor.contig_start - candidate.contig_end;
            } 
        }
        
        if (static_cast<int>(contig_distance) < min_alignment_distance) {
            return true;
        }
    }
    
    return false;
}

void Segmentation::process_read(uint32_t read_idx,
                                const vector<size_t>& alignment_indices,
                                const vector<Alignment>& alignments,
                                vector<ReadBreakpoint>& output_breakpoints)
{
    string read_id = store.get_read_id(read_idx);
    uint32_t read_length = store.get_reads()[read_idx].length;
    
    for (size_t i = 0; i < alignment_indices.size(); ++i) {
        const Alignment& anchor = alignments[alignment_indices[i]];
        
        if (!is_valid_anchor(anchor)) {
            continue;
        }
        
        // check both ends of the anchor
        // for positive strand: end -> next, start -> previous
        // for negative strand: end -> previous, start -> next
        for (bool check_end : {true, false}) {
            bool check_next = (anchor.is_reverse && !check_end) || (!anchor.is_reverse && check_end);
            
            if (is_too_close_to_read_boundary(anchor, read_length, check_next)) {
                continue;
            }
            
            if (has_nearby_same_contig_alignment(anchor, i, check_next, alignment_indices, alignments)) {
                continue;
            }
            
            string breakpoint_type = get_breakpoint_type(anchor.is_reverse, check_end);
            uint32_t coord = get_breakpoint_coord(anchor, check_end);
            string contig_id = store.get_contig_id(anchor.contig_index);
            
            uint32_t anchor_length = anchor.read_end - anchor.read_start;
            uint32_t anchor_mutations = anchor.get_mutation_count();
            
            ReadBreakpoint bp(current_lib_id, read_id, contig_id, coord, breakpoint_type,
                            anchor_length, anchor_mutations, 0);
            output_breakpoints.push_back(bp);
        }
    }
}

void Segmentation::detect_breakpoints(const string& lib_id, 
                                      vector<ReadBreakpoint>& output_breakpoints)
{
    cout << "detecting breakpoints for library " << lib_id << "..." << endl;
    current_lib_id = lib_id;
    
    output_breakpoints.clear();
    
    const auto& reads = store.get_reads();
    
    // ensure read-to-alignment index is built
    massert(store.is_read_alignment_index_built(), 
           "read-to-alignment index must be built before detecting breakpoints");
    
    const auto& read_to_alignments = store.get_read_to_alignment_indices();
    const auto& alignments = store.get_alignments();
    
    size_t tested_reads = 0;
    size_t reads_with_breakpoints = 0;
    
    // process each read
    for (size_t read_idx = 0; read_idx < reads.size(); ++read_idx) {
        auto it = read_to_alignments.find(static_cast<uint32_t>(read_idx));
        if (it == read_to_alignments.end()) {
            continue;
        }
        
        const vector<size_t>& alignment_indices = it->second;
        
        size_t before_count = output_breakpoints.size();
        process_read(static_cast<uint32_t>(read_idx), alignment_indices, alignments, output_breakpoints);
        
        if (output_breakpoints.size() > before_count) {
            reads_with_breakpoints++;
        }
        
        tested_reads++;
    }
    
    cout << "tested " << tested_reads << " reads" << endl;
    cout << "found " << output_breakpoints.size() << " breakpoints in " 
         << reads_with_breakpoints << " reads" << endl;
}

