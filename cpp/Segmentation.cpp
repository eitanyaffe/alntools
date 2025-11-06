#include "Segmentation.h"
#include <algorithm>
#include <iostream>
#include <set>

using namespace std;

Segmentation::Segmentation(const AlignmentStore& store,
                         int /* max_margin */,
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
        
        if (alignment_indices.size() < 2) {
            continue;
        }
        
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

bool Segmentation::is_valid_dangle(const Alignment& dangle, const Alignment& anchor) const
{
    uint32_t dangle_length = dangle.read_end - dangle.read_start;
    if (static_cast<int>(dangle_length) < min_dangle_length) {
        return false;
    }
    
    if (!is_valid_dangle_position(anchor, dangle)) {
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

bool Segmentation::is_valid_dangle_position(const Alignment& anchor, const Alignment& dangle) const
{
    // different contig - valid dangle
    if (anchor.contig_index != dangle.contig_index) {
        return true;
    }
    
    // same contig - check distance
    uint32_t distance;
    if (anchor.contig_end <= dangle.contig_start) {
        distance = dangle.contig_start - anchor.contig_end;
    } else if (dangle.contig_end <= anchor.contig_start) {
        distance = anchor.contig_start - dangle.contig_end;
    } else {
        // overlapping - not valid
        return false;
    }
    
    return static_cast<int>(distance) >= min_alignment_distance;
}

void Segmentation::process_read(uint32_t read_idx,
                                const vector<size_t>& alignment_indices,
                                const vector<Alignment>& alignments,
                                vector<ReadBreakpoint>& output_breakpoints)
{
    string read_id = store.get_read_id(read_idx);
    
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
            int dangle_idx = check_next ? (i + 1) : (i - 1);
            
            if (dangle_idx < 0 || dangle_idx >= static_cast<int>(alignment_indices.size())) {
                continue;
            }
            
            const Alignment& dangle = alignments[alignment_indices[dangle_idx]];
            if (!is_valid_dangle(dangle, anchor)) {
                continue;
            }
            
            string breakpoint_type = get_breakpoint_type(anchor.is_reverse, check_end);
            uint32_t coord = get_breakpoint_coord(anchor, check_end);
            string contig_id = store.get_contig_id(anchor.contig_index);
            
            uint32_t anchor_length = anchor.read_end - anchor.read_start;
            uint32_t anchor_mutations = anchor.get_mutation_count();
            uint32_t dangle_length = dangle.read_end - dangle.read_start;
            
            ReadBreakpoint bp(current_lib_id, read_id, contig_id, coord, breakpoint_type,
                            anchor_length, anchor_mutations, dangle_length);
            output_breakpoints.push_back(bp);
        }
    }
}

