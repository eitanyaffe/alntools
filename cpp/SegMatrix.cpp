#include "SegMatrix.h"
#include "utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////
// constructor
////////////////////////////////////////////////////////////////////////////////////////////////////

SegMatrix::SegMatrix(AlignmentStore& store_ref) : store(store_ref), has_clusters(false), debug_read_count(0) {}

////////////////////////////////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////////////////////////////////

void SegMatrix::load_segments(const string& ifn_segments)
{
    ifstream file(ifn_segments);
    massert(file.is_open(), "could not open file %s", ifn_segments.c_str());
    
    string line;
    int seg_ind = -1, contig_ind = -1, start_ind = -1, end_ind = -1;
    
    segment_index_map.clear();
    segments.clear();

    // parse header
    if (getline(file, line)) {
        istringstream header_iss(line);
        string header;
        int col_index = 0;
        while (getline(header_iss, header, '\t')) {
            if (header == "segment") seg_ind = col_index;
            else if (header == "contig") contig_ind = col_index;
            else if (header == "start") start_ind = col_index;
            else if (header == "end") end_ind = col_index;
            col_index++;
        }
        
        massert(seg_ind >= 0, "column 'segment' not found in %s", ifn_segments.c_str());
        massert(contig_ind >= 0, "column 'contig' not found in %s", ifn_segments.c_str());
        massert(start_ind >= 0, "column 'start' not found in %s", ifn_segments.c_str());
        massert(end_ind >= 0, "column 'end' not found in %s", ifn_segments.c_str());
    }
    
    // parse data rows - keep all segments (including short ones)
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        istringstream iss(line);
        vector<string> fields;
        string field;
        while (getline(iss, field, '\t')) {
            fields.push_back(field);
        }
        
        if (fields.size() == 0) continue;
            
        SegMatrixSegment seg;
        seg.id = fields[seg_ind];
        seg.contig = fields[contig_ind];
        seg.start = atoi(fields[start_ind].c_str());
        seg.end = atoi(fields[end_ind].c_str());
        seg.length = seg.end - seg.start + 1;
        seg.index = segments.size();

        uint32_t seg_start_0based = seg.start - 1;
        uint32_t seg_end_0based = seg.end;
        if (seg.length >= min_segment_length) {
            uint32_t left_start = seg_start_0based + side_margin;
            uint32_t left_end = left_start + side_length;
            uint32_t right_end = seg_end_0based > side_margin ? seg_end_0based - side_margin : seg_end_0based;
            uint32_t right_start = right_end > side_length ? right_end - side_length : seg_start_0based;
            seg.left_side_start = left_start;
            seg.left_side_end = left_end;
            seg.right_side_start = right_start;
            seg.right_side_end = right_end;
        } else {
            seg.left_side_start = seg.left_side_end = seg_start_0based;
            seg.right_side_start = seg.right_side_end = seg_end_0based;
        }

        segments.push_back(seg);
        segment_index_map[seg.id] = seg.index;
    }
}

void SegMatrix::load_cluster_mapping(const string& ifn_clusters)
{
    if (ifn_clusters.empty()) {
        has_clusters = false;
        return;
    }
    
    ifstream file(ifn_clusters);
    massert(file.is_open(), "could not open file %s", ifn_clusters.c_str());
    
    string line;
    int seg_ind = -1, cluster_ind = -1, strand_ind = -1;
    
    segment_to_cluster.clear();
    segment_to_strand.clear();
    
    // parse header
    if (getline(file, line)) {
        istringstream header_iss(line);
        string header;
        int col_index = 0;
        while (getline(header_iss, header, '\t')) {
            if (header == "segment") seg_ind = col_index;
            else if (header == "csegment") cluster_ind = col_index;
            else if (header == "strand") strand_ind = col_index;
            col_index++;
        }
        
        massert(seg_ind >= 0, "column 'segment' not found in %s", ifn_clusters.c_str());
        massert(cluster_ind >= 0, "column 'csegment' not found in %s", ifn_clusters.c_str());
        massert(strand_ind >= 0, "column 'strand' not found in %s", ifn_clusters.c_str());
    }
    
    // parse data rows
    while (getline(file, line)) {
        if (line.empty()) continue;
        
        istringstream iss(line);
        vector<string> fields;
        string field;
        while (getline(iss, field, '\t')) {
            fields.push_back(field);
        }
        
        if (fields.size() == 0) continue;
        
        string seg_id = fields[seg_ind];
        string cluster_id = fields[cluster_ind];
        string strand = fields[strand_ind];
        
        segment_to_cluster[seg_id] = cluster_id;
        segment_to_strand[seg_id] = strand;
    }
    
    has_clusters = true;
    cout << "loaded cluster mapping for " << segment_to_cluster.size() << " segments" << endl;
}

void SegMatrix::load_library()
{
    store.count_short_indels(min_indel_length);
    
    if (!store.is_read_alignment_index_built()) {
        store.init_read_alignment_index();
    }
    
    cout << "using alignment store with " << store.get_read_count() << " reads and " 
         << store.get_alignment_count() << " alignments" << endl;
}

void SegMatrix::build_contig_to_segments_map()
{
    contig_to_segments.clear();
    for (size_t i = 0; i < segments.size(); ++i) {
        contig_to_segments[segments[i].contig].push_back(i);
    }
    cout << "built contig to segments map for " << contig_to_segments.size() << " contigs" << endl;
}

void SegMatrix::write_matrix(const map<tuple<string, string, string, string>, uint32_t>& matrix,
                              const string& filename,
                              const string& matrix_name) const
{
    ofstream out(filename);
    massert(out.is_open(), "could not open file %s", filename.c_str());
    
    out << "seg_src\tseg_tgt\tside_src\tside_tgt\tcontig_src\tstart_src\tend_src\tcontig_tgt\tstart_tgt\tend_tgt\ttotal_read_count\tassociated_read_count\tcount" << endl;
    
    for (const auto& entry : matrix) {
        const auto& key = entry.first;
        uint32_t count = entry.second;
        const string& seg_src = get<0>(key);
        const string& seg_tgt = get<1>(key);
        const string& side_src = get<2>(key);
        const string& side_tgt = get<3>(key);
        auto idx_src_it = segment_index_map.find(seg_src);
        auto idx_tgt_it = segment_index_map.find(seg_tgt);
        massert(idx_src_it != segment_index_map.end(), "segment %s not found", seg_src.c_str());
        massert(idx_tgt_it != segment_index_map.end(), "segment %s not found", seg_tgt.c_str());
        const SegMatrixSegment& seg_info_src = segments[idx_src_it->second];
        const SegMatrixSegment& seg_info_tgt = segments[idx_tgt_it->second];
        
        uint64_t total_src = 0;
        uint64_t assoc_src = 0;
        auto stats_src_it = segment_stats.find(seg_src);
        if (stats_src_it != segment_stats.end()) {
            const SegmentStats& stats = stats_src_it->second;
            total_src = (side_src == "left") ? stats.total_read_count_left : stats.total_read_count_right;
            assoc_src = (side_src == "left") ? stats.associated_read_count_left : stats.associated_read_count_right;
        }
        out << seg_src << "\t" << seg_tgt << "\t"
            << side_src << "\t" << side_tgt << "\t"
            << seg_info_src.contig << "\t" << seg_info_src.start << "\t" << seg_info_src.end << "\t"
            << seg_info_tgt.contig << "\t" << seg_info_tgt.start << "\t" << seg_info_tgt.end << "\t"
            << total_src << "\t" << assoc_src << "\t"
            << count << endl;
        
        if (seg_src != seg_tgt) {
            uint64_t total_tgt = 0;
            uint64_t assoc_tgt = 0;
            auto stats_tgt_it = segment_stats.find(seg_tgt);
            if (stats_tgt_it != segment_stats.end()) {
                const SegmentStats& stats = stats_tgt_it->second;
                total_tgt = (side_tgt == "left") ? stats.total_read_count_left : stats.total_read_count_right;
                assoc_tgt = (side_tgt == "left") ? stats.associated_read_count_left : stats.associated_read_count_right;
            }
            out << seg_tgt << "\t" << seg_src << "\t"
                << side_tgt << "\t" << side_src << "\t"
                << seg_info_tgt.contig << "\t" << seg_info_tgt.start << "\t" << seg_info_tgt.end << "\t"
                << seg_info_src.contig << "\t" << seg_info_src.start << "\t" << seg_info_src.end << "\t"
                << total_tgt << "\t" << assoc_tgt << "\t"
                << count << endl;
        }
    }
    
    out.close();
    cout << "wrote " << matrix_name << ": " << filename << " (" 
         << matrix.size() << " entries)" << endl;
}

void SegMatrix::write_matrices(const string& odir)
{
    write_matrix(adjacency_matrix, odir + "/adjacency_matrix.txt", "adjacency matrix");
    write_matrix(reach_matrix, odir + "/reach_matrix.txt", "reach matrix");
}

void SegMatrix::write_segment_summary(const string& odir)
{
    string summary_file = odir + "/segment_summary.txt";
    
    ofstream out(summary_file);
    massert(out.is_open(), "could not open file %s", summary_file.c_str());
    
    out << "segment\tcontig\tstart\tend\tlength\tis_short\ttotal_read_count_left\ttotal_read_count_right\tassociated_read_count_left\tassociated_read_count_right\tassociated_segment_count_left\tassociated_segment_count_right" << endl;
    
    // iterate through all segments to ensure we output all of them
    for (const auto& seg : segments) {
        const string& seg_id = seg.id;
        auto stats_it = segment_stats.find(seg_id);
        
        uint64_t total_left = 0;
        uint64_t total_right = 0;
        uint64_t assoc_left = 0;
        uint64_t assoc_right = 0;
        uint64_t assoc_seg_left = 0;
        uint64_t assoc_seg_right = 0;
        
        if (stats_it != segment_stats.end()) {
            const SegmentStats& stats = stats_it->second;
            total_left = stats.total_read_count_left;
            total_right = stats.total_read_count_right;
            assoc_left = stats.associated_read_count_left;
            assoc_right = stats.associated_read_count_right;
            assoc_seg_left = stats.associated_segments_left.size();
            assoc_seg_right = stats.associated_segments_right.size();
        }
        
        bool is_short = (seg.length < min_segment_length);
        
        out << seg_id << "\t" << seg.contig << "\t" << seg.start << "\t" << seg.end 
            << "\t" << seg.length << "\t" << (is_short ? "T" : "F")
            << "\t" << total_left << "\t" << total_right
            << "\t" << assoc_left << "\t" << assoc_right
            << "\t" << assoc_seg_left << "\t" << assoc_seg_right
            << endl;
    }
    
    out.close();
    cout << "wrote segment summary: " << summary_file << " (" 
         << segments.size() << " segments)" << endl;
}

void SegMatrix::write_read_details(const string& odir)
{
    string adjacency_file = odir + "/read_associations_adjacency.txt";
    ofstream out_adj(adjacency_file);
    massert(out_adj.is_open(), "could not open file %s", adjacency_file.c_str());
    
    out_adj << "seg_src\tside_src\tseg_tgt\tside_tgt\tread_id\t"
            << "aln_src_contig\taln_src_read_start\taln_src_read_end\taln_src_contig_start\taln_src_contig_end\taln_src_is_reverse\t"
            << "aln_tgt_contig\taln_tgt_read_start\taln_tgt_read_end\taln_tgt_contig_start\taln_tgt_contig_end\taln_tgt_is_reverse" << endl;
    for (const auto& entry : adjacency_reads) {
        out_adj << entry.seg_src << "\t" << entry.side_src << "\t"
                << entry.seg_tgt << "\t" << entry.side_tgt << "\t"
                << entry.read_id << "\t"
                << entry.aln_src_contig << "\t"
                << (entry.aln_src_read_start + 1) << "\t" << entry.aln_src_read_end << "\t"
                << (entry.aln_src_contig_start + 1) << "\t" << entry.aln_src_contig_end << "\t"
                << (entry.aln_src_is_reverse ? "-" : "+") << "\t"
                << entry.aln_tgt_contig << "\t"
                << (entry.aln_tgt_read_start + 1) << "\t" << entry.aln_tgt_read_end << "\t"
                << (entry.aln_tgt_contig_start + 1) << "\t" << entry.aln_tgt_contig_end << "\t"
                << (entry.aln_tgt_is_reverse ? "-" : "+") << endl;
    }
    out_adj.close();
    cout << "wrote adjacency read details: " << adjacency_file << " (" 
         << adjacency_reads.size() << " entries)" << endl;
    
    string reach_file = odir + "/read_associations_reach.txt";
    ofstream out_reach(reach_file);
    massert(out_reach.is_open(), "could not open file %s", reach_file.c_str());
    
    out_reach << "seg_src\tside_src\tseg_tgt\tside_tgt\tread_id\t"
              << "aln_src_contig\taln_src_read_start\taln_src_read_end\taln_src_contig_start\taln_src_contig_end\taln_src_is_reverse\t"
              << "aln_tgt_contig\taln_tgt_read_start\taln_tgt_read_end\taln_tgt_contig_start\taln_tgt_contig_end\taln_tgt_is_reverse" << endl;
    for (const auto& entry : reach_reads) {
        out_reach << entry.seg_src << "\t" << entry.side_src << "\t"
                  << entry.seg_tgt << "\t" << entry.side_tgt << "\t"
                  << entry.read_id << "\t"
                  << entry.aln_src_contig << "\t"
                  << (entry.aln_src_read_start + 1) << "\t" << entry.aln_src_read_end << "\t"
                  << (entry.aln_src_contig_start + 1) << "\t" << entry.aln_src_contig_end << "\t"
                  << (entry.aln_src_is_reverse ? "-" : "+") << "\t"
                  << entry.aln_tgt_contig << "\t"
                  << (entry.aln_tgt_read_start + 1) << "\t" << entry.aln_tgt_read_end << "\t"
                  << (entry.aln_tgt_contig_start + 1) << "\t" << entry.aln_tgt_contig_end << "\t"
                  << (entry.aln_tgt_is_reverse ? "-" : "+") << endl;
    }
    out_reach.close();
    cout << "wrote reach read details: " << reach_file << " (" 
         << reach_reads.size() << " entries)" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions called by process_read
////////////////////////////////////////////////////////////////////////////////////////////////////

vector<const Alignment*> SegMatrix::get_sorted_alignments_for_read(uint32_t read_index)
{
    vector<const Alignment*> result;
    
    const auto& read_to_alignment_indices = store.get_read_to_alignment_indices();
    auto it = read_to_alignment_indices.find(read_index);
    if (it == read_to_alignment_indices.end()) {
        return result;
    }
    
    const vector<size_t>& alignment_indices = it->second;
    const vector<Alignment>& all_alignments = store.get_alignments();
    
    for (size_t idx : alignment_indices) {
        result.push_back(&all_alignments[idx]);
    }
    
    return result;
}

vector<const Alignment*> SegMatrix::filter_overlapping_alignments(const vector<const Alignment*>& alignments, bool debug_mode)
{
    if (alignments.empty()) {
        return alignments;
    }
    
    vector<const Alignment*> filtered;
    filtered.push_back(alignments[0]);
    
    for (size_t i = 1; i < alignments.size(); ++i) {
        const Alignment* prev = filtered.back();
        const Alignment* curr = alignments[i];
        
        // calculate overlap in read coordinates
        uint32_t overlap_start = max(prev->read_start, curr->read_start);
        uint32_t overlap_end = min(prev->read_end, curr->read_end);
        
        if (overlap_end > overlap_start) {
            uint32_t overlap = overlap_end - overlap_start;
            if (overlap > max_margin) {
                if (debug_mode) {
                    cout << "  alignment[" << i << "] was filtered because overlap=" << overlap 
                         << " > max_margin=" << max_margin << endl;
                }
                continue; // skip this alignment
            }
        }
        
        filtered.push_back(curr);
    }
    
    return filtered;
}

bool SegMatrix::is_better_alignment(const Alignment* a, const Alignment* b)
{
    if (!a) return false;
    if (!b) return true;
    
    uint32_t a_read_span = a->read_end - a->read_start;
    uint32_t b_read_span = b->read_end - b->read_start;
    if (a_read_span != b_read_span) {
        return a_read_span > b_read_span;
    }
    
    uint32_t a_contig_span = a->contig_end - a->contig_start;
    uint32_t b_contig_span = b->contig_end - b->contig_start;
    if (a_contig_span != b_contig_span) {
        return a_contig_span > b_contig_span;
    }
    
    if (a->mutations.size() != b->mutations.size()) {
        return a->mutations.size() < b->mutations.size();
    }
    
    if (a->read_start != b->read_start) {
        return a->read_start < b->read_start;
    }
    
    return a->contig_index < b->contig_index;
}

vector<const Alignment*> SegMatrix::filter_candidates_by_quality(const vector<const Alignment*>& alignments)
{
    vector<const Alignment*> filtered;
    
    for (const Alignment* aln : alignments) {
        uint32_t contig_span = aln->contig_end - aln->contig_start;
        if (contig_span < side_length) {
            continue;
        }
        
        double mutation_percent = compute_mutation_percent(*aln, aln->contig_start, aln->contig_end);
        if (mutation_percent > max_mutation_percent) {
            continue;
        }
        
        filtered.push_back(aln);
    }
    
    return filtered;
}

bool SegMatrix::has_excessive_overlap(const Alignment* aln, const vector<const Alignment*>& selected)
{
    for (const Alignment* sel : selected) {
        uint32_t overlap_start = max(aln->read_start, sel->read_start);
        uint32_t overlap_end = min(aln->read_end, sel->read_end);
        
        if (overlap_end > overlap_start) {
            uint32_t overlap = overlap_end - overlap_start;
            if (overlap > max_margin) {
                return true;
            }
        }
    }
    return false;
}

const Alignment* SegMatrix::find_top_alignment(const vector<const Alignment*>& candidates,
                                               const map<uint32_t, uint32_t>& contig_to_read_length)
{
    if (candidates.empty()) {
        return nullptr;
    }
    
    vector<pair<uint32_t, uint32_t>> sorted_contigs;
    for (const auto& kv : contig_to_read_length) {
        sorted_contigs.push_back({kv.first, kv.second});
    }
    sort(sorted_contigs.begin(), sorted_contigs.end(),
         [](const pair<uint32_t, uint32_t>& a, const pair<uint32_t, uint32_t>& b) {
             return a.second > b.second;
         });
    
    for (const auto& contig_pair : sorted_contigs) {
        uint32_t contig_index = contig_pair.first;
        const Alignment* best = nullptr;
        
        for (const Alignment* cand : candidates) {
            if (cand->contig_index == contig_index) {
                if (is_better_alignment(cand, best)) {
                    best = cand;
                }
            }
        }
        
        if (best) {
            return best;
        }
    }
    
    const Alignment* best = nullptr;
    for (const Alignment* cand : candidates) {
        if (is_better_alignment(cand, best)) {
            best = cand;
        }
    }
    
    return best;
}

vector<const Alignment*> SegMatrix::get_parsimony_read_alignments(const vector<const Alignment*>& alignments)
{
    if (alignments.empty()) {
        return alignments;
    }
    
    vector<const Alignment*> candidates = filter_candidates_by_quality(alignments);
    vector<const Alignment*> selected;
    map<uint32_t, uint32_t> contig_to_read_length;
    
    while (!candidates.empty()) {
        const Alignment* top = find_top_alignment(candidates, contig_to_read_length);
        
        if (!top) {
            break;
        }
        
        if (!has_excessive_overlap(top, selected)) {
            selected.push_back(top);
            uint32_t read_span = top->read_end - top->read_start;
            contig_to_read_length[top->contig_index] += read_span;
        }
        
        auto it = find(candidates.begin(), candidates.end(), top);
        if (it != candidates.end()) {
            candidates.erase(it);
        }
    }
    
    return selected;
}

vector<ReadInterval> SegMatrix::intersect_alignment_with_segments(const Alignment& aln)
{
    vector<ReadInterval> intervals;
    
    string contig_id = store.get_contig_id(aln.contig_index);
    auto contig_segments_it = contig_to_segments.find(contig_id);
    if (contig_segments_it == contig_to_segments.end()) {
        return intervals;
    }
    
    const vector<size_t>& seg_indices = contig_segments_it->second;
    
    // convert segment coordinates to 0-based half-open for intersection
    uint32_t aln_contig_start = aln.contig_start;
    uint32_t aln_contig_end = aln.contig_end;
    
    for (size_t seg_idx : seg_indices) {
        const SegMatrixSegment& seg = segments[seg_idx];
        
        // convert segment to 0-based half-open
        uint32_t seg_start_0based = seg.start - 1;
        uint32_t seg_end_0based = seg.end;
        
        // check intersection
        uint32_t intersection_start = max(aln_contig_start, seg_start_0based);
        uint32_t intersection_end = min(aln_contig_end, seg_end_0based);
        
        if (intersection_end > intersection_start) {
            ReadInterval interval;
            interval.segment_id = seg.id;
            interval.contig = contig_id;
            interval.contig_start = intersection_start;
            interval.contig_end = intersection_end;
            interval.is_reverse = aln.is_reverse;
            interval.is_short_segment = (seg.length < min_segment_length);
            interval.segment_index = seg_idx;
            
            // use precise coordinate mapping by processing mutations
            uint32_t read_start_precise = contig_to_read_coord(aln, intersection_start, store);
            uint32_t read_end_precise = contig_to_read_coord(aln, intersection_end, store);
            
            // ensure proper ordering: start is min, end is max
            interval.read_start_est = std::min(read_start_precise, read_end_precise);
            interval.read_end_est = std::max(read_start_precise, read_end_precise);
            
            intervals.push_back(interval);
        }
    }

    if (aln.is_reverse) {
        std::reverse(intervals.begin(), intervals.end());
    }

    return intervals;
}

double SegMatrix::compute_mutation_percent(const Alignment& aln, uint32_t contig_start, uint32_t contig_end)
{
    uint32_t mutation_count = 0;
    uint32_t region_length = contig_end - contig_start;
    
    if (region_length == 0) {
        return 0.0;
    }
    
    for (uint32_t mutation_index : aln.mutations) {
        if (should_skip_short_indel(store, aln.contig_index, mutation_index)) {
            continue;
        }
        
        const Mutation& mutation = store.get_mutation(aln.contig_index, mutation_index);
        uint32_t mut_pos = mutation.position;
        
        if (mut_pos >= contig_start && mut_pos < contig_end) {
            mutation_count++;
        }
    }
    
    return (static_cast<double>(mutation_count) / static_cast<double>(region_length)) * 100.0;
}

void SegMatrix::compute_side_coverage_and_mutations(ReadInterval& interval, const Alignment& aln)
{
    // get segment using stored index
    if (interval.segment_index >= segments.size()) {
        return;
    }

    const SegMatrixSegment& seg = segments[interval.segment_index];

    if (seg.length < min_segment_length) {
        return;
    }

    // convert segment to 0-based half-open
    uint32_t seg_start_0based = seg.start - 1;
    uint32_t seg_end_0based = seg.end;

    uint32_t left_side_start = seg.left_side_start;
    uint32_t left_side_end = seg.left_side_end;
    uint32_t right_side_start = seg.right_side_start;
    uint32_t right_side_end = seg.right_side_end;

    bool left_region_valid = left_side_end > left_side_start && left_side_end <= seg_end_0based;
    bool right_region_valid = right_side_end > right_side_start && right_side_start >= seg_start_0based;

    if (left_region_valid && interval.contig_start <= left_side_start && interval.contig_end >= left_side_end) {
        interval.covers_left_side = true;
        interval.left_mutation_percent = compute_mutation_percent(aln, left_side_start, left_side_end);
    }

    if (right_region_valid && interval.contig_start <= right_side_start && interval.contig_end >= right_side_end) {
        interval.covers_right_side = true;
        interval.right_mutation_percent = compute_mutation_percent(aln, right_side_start, right_side_end);
    }
}

bool SegMatrix::is_pointing_out(const ReadInterval& interval) const
{
    if (interval.is_reverse) {
        return interval.covers_left_side;
    } else {
        return interval.covers_right_side;
    }
}

bool SegMatrix::is_pointing_in(const ReadInterval& interval) const
{
    if (interval.is_reverse) {
        return interval.covers_right_side;
    } else {
        return interval.covers_left_side;
    }
}

bool SegMatrix::is_within_adjacency_distance(const ReadInterval& interval1, const ReadInterval& interval2) const
{
    // intervals are adjacent if the gap between their estimated read coordinates is small
    uint32_t start1 = interval1.read_start_est;
    uint32_t end1 = interval1.read_end_est;
    uint32_t start2 = interval2.read_start_est;
    uint32_t end2 = interval2.read_end_est;

    if (start2 < start1) {
        std::swap(start1, start2);
        std::swap(end1, end2);
    }

    uint64_t gap = 0;
    if (end1 < start2) {
        gap = static_cast<uint64_t>(start2) - static_cast<uint64_t>(end1);
    } else if (end2 < start1) {
        gap = static_cast<uint64_t>(start1) - static_cast<uint64_t>(end2);
    }

    return gap <= static_cast<uint64_t>(max_adjacency_distance);
}

void SegMatrix::print_read_debug_info(uint32_t read_index, const string& read_id,
                                      const vector<const Alignment*>& alignments_before,
                                      const vector<const Alignment*>& alignments_after,
                                      const vector<ReadInterval>& intervals) const
{
    cout << "================================================================================" << endl;
    cout << "read[" << read_index << "] " << read_id << ": " << alignments_before.size() << " alignments before filtering" << endl;
    cout << "================================================================================" << endl;
    for (size_t i = 0; i < alignments_before.size(); ++i) {
        const Alignment* aln = alignments_before[i];
        string contig_id = store.get_contig_id(aln->contig_index);
        cout << "  aln[" << i << "]: contig=" << contig_id 
             << ", read_coords=[" << aln->read_start << "-" << aln->read_end << ")"
             << ", contig_coords=[" << aln->contig_start << "-" << aln->contig_end << ")"
             << ", reverse=" << (aln->is_reverse ? "T" : "F") << endl;
    }
    
    cout << "read[" << read_index << "] " << read_id << ": " << alignments_after.size() << " alignments after filtering" << endl;
    for (size_t i = 0; i < alignments_after.size(); ++i) {
        const Alignment* aln = alignments_after[i];
        string contig_id = store.get_contig_id(aln->contig_index);
        cout << "  aln[" << i << "]: contig=" << contig_id 
             << ", read_coords=[" << aln->read_start << "-" << aln->read_end << ")"
             << ", contig_coords=[" << aln->contig_start << "-" << aln->contig_end << ")"
             << ", reverse=" << (aln->is_reverse ? "T" : "F") << endl;
    }
    
    print_interval_sequence(intervals);
}

void SegMatrix::print_read_header(uint32_t read_index, const string& read_id) const
{
    cout << "\n=== DEBUG READ #" << (debug_read_count + 1) << " ===" << endl;
    cout << "read_index: " << read_index << ", read_id: " << read_id << endl;
}

void SegMatrix::print_interval(const ReadInterval& interval, size_t idx) const
{
    // get segment info for debug output
    massert(interval.segment_index < segments.size(),
           "interval[%zu] has invalid segment_index %zu (segments.size()=%zu) for segment_id=%s",
           idx, interval.segment_index, segments.size(), interval.segment_id.c_str());
    
    const SegMatrixSegment& seg = segments[interval.segment_index];
    
    uint32_t interval_length = interval.contig_end - interval.contig_start;
    cout << "  interval[" << idx << "]: segment=" << interval.segment_id
         << " contig=" << interval.contig
         << " seg_coords=[" << (seg.start - 1) << "-" << seg.end << "]"
         << " seg_len=" << seg.length
         << ", interval_contig=[" << interval.contig_start << "-" << interval.contig_end << ")"
         << " interval_len=" << interval_length
         << ", read_est=[" << interval.read_start_est << "-" << interval.read_end_est << ")"
         << ", reverse=" << (interval.is_reverse ? "T" : "F")
         << ", left=" << (interval.covers_left_side ? "T" : "F")
         << ", right=" << (interval.covers_right_side ? "T" : "F")
         << ", left_mut%=" << interval.left_mutation_percent
         << ", right_mut%=" << interval.right_mutation_percent
         << ", short=" << (interval.is_short_segment ? "T" : "F") << endl;
}

void SegMatrix::print_interval_sequence(const vector<ReadInterval>& intervals) const
{
    cout << "intervals (" << intervals.size() << " total):" << endl;
    for (size_t i = 0; i < intervals.size(); ++i) {
        print_interval(intervals[i], i);
    }
}

void SegMatrix::print_short_segment_skipped(const ReadInterval& interval) const
{
    cout << "  -> skipped short segment: " << interval.segment_id << endl;
}


void SegMatrix::print_exit_detected(const ReadInterval& interval, double mutation_percent) const
{
    string side = interval.is_reverse ? "left" : "right";
    cout << "  -> EXIT detected: leaving segment " << interval.segment_id 
         << " on " << side << " side (mutation%=" << mutation_percent << ")" << endl;
}

void SegMatrix::print_entry_detected(const ReadInterval& interval, double mutation_percent) const
{
    string side = interval.is_reverse ? "right" : "left";
    cout << "  -> ENTRY detected: entering segment " << interval.segment_id 
         << " on " << side << " side (mutation%=" << mutation_percent << ")" << endl;
}

void SegMatrix::print_association(const ReadInterval& exit_interval, const ReadInterval& entry_interval, 
                                  bool is_adjacent) const
{
    string exit_side = exit_interval.is_reverse ? "left" : "right";
    string entry_side = entry_interval.is_reverse ? "right" : "left";
    cout << "  -> ASSOCIATION: " << exit_interval.segment_id << "(" << exit_side << ") -> " 
         << entry_interval.segment_id << "(" << entry_side << ")";
    if (is_adjacent) {
        cout << " [ADJACENT]";
    }
    cout << endl;
}

bool SegMatrix::should_skip_mutation_threshold(double mutation_percent, const ReadInterval& interval,
                                                const string& direction, bool debug_mode) const
{
    if (mutation_percent > max_mutation_percent) {
        if (debug_mode) {
            string side;
            string preposition;
            if (direction == "exit") {
                side = interval.is_reverse ? "left" : "right";
                preposition = "from";
            } else {
                side = interval.is_reverse ? "right" : "left";
                preposition = "to";
            }
            cout << "  -> skipping " << direction << " " << preposition << " " << interval.segment_id 
                 << " (" << side << " side): mutation%=" << mutation_percent 
                 << " > threshold=" << max_mutation_percent << endl;
        }
        return true;
    }
    return false;
}

void SegMatrix::create_association(const ReadInterval& exit_interval, const ReadInterval& entry_interval,
                                   bool debug_mode,
                                   std::unordered_map<std::string, SegmentReadContribution>& read_contrib,
                                   const std::string& read_id)
{
    string exit_side = exit_interval.is_reverse ? "left" : "right";
    string entry_side = entry_interval.is_reverse ? "right" : "left";

    SegmentStats& exit_stats = segment_stats[exit_interval.segment_id];
    SegmentStats& entry_stats = segment_stats[entry_interval.segment_id];

    if (exit_side == "left") {
        exit_stats.associated_segments_left.insert(entry_interval.segment_id);
    } else {
        exit_stats.associated_segments_right.insert(entry_interval.segment_id);
    }
    if (entry_side == "left") {
        entry_stats.associated_segments_left.insert(exit_interval.segment_id);
    } else {
        entry_stats.associated_segments_right.insert(exit_interval.segment_id);
    }

    SegmentReadContribution& exit_contrib = read_contrib[exit_interval.segment_id];
    if (exit_side == "left") {
        exit_contrib.associated_left = true;
    } else {
        exit_contrib.associated_right = true;
    }
    SegmentReadContribution& entry_contrib = read_contrib[entry_interval.segment_id];
    if (entry_side == "left") {
        entry_contrib.associated_left = true;
    } else {
        entry_contrib.associated_right = true;
    }

    std::string seg_src = exit_interval.segment_id;
    std::string seg_tgt = entry_interval.segment_id;
    std::string side_src = exit_side;
    std::string side_tgt = entry_side;

    auto key = canonize_key(seg_src, seg_tgt, side_src, side_tgt);

    bool is_adjacent = is_within_adjacency_distance(exit_interval, entry_interval);
    if (debug_mode) {
        print_association(exit_interval, entry_interval, is_adjacent);
    }

    reach_matrix[key]++;
    if (output_read_details_flag) {
        ReadAssociation assoc;
        assoc.seg_src = seg_src;
        assoc.side_src = side_src;
        assoc.seg_tgt = seg_tgt;
        assoc.side_tgt = side_tgt;
        assoc.read_id = read_id;
        assoc.aln_src_contig = exit_interval.contig;
        assoc.aln_src_read_start = exit_interval.read_start_est;
        assoc.aln_src_read_end = exit_interval.read_end_est;
        assoc.aln_src_contig_start = exit_interval.contig_start;
        assoc.aln_src_contig_end = exit_interval.contig_end;
        assoc.aln_src_is_reverse = exit_interval.is_reverse;
        assoc.aln_tgt_contig = entry_interval.contig;
        assoc.aln_tgt_read_start = entry_interval.read_start_est;
        assoc.aln_tgt_read_end = entry_interval.read_end_est;
        assoc.aln_tgt_contig_start = entry_interval.contig_start;
        assoc.aln_tgt_contig_end = entry_interval.contig_end;
        assoc.aln_tgt_is_reverse = entry_interval.is_reverse;
        reach_reads.push_back(assoc);
    }

    if (is_adjacent) {
        adjacency_matrix[key]++;
        if (output_read_details_flag) {
            ReadAssociation assoc;
            assoc.seg_src = seg_src;
            assoc.side_src = side_src;
            assoc.seg_tgt = seg_tgt;
            assoc.side_tgt = side_tgt;
            assoc.read_id = read_id;
            assoc.aln_src_contig = exit_interval.contig;
            assoc.aln_src_read_start = exit_interval.read_start_est;
            assoc.aln_src_read_end = exit_interval.read_end_est;
            assoc.aln_src_contig_start = exit_interval.contig_start;
            assoc.aln_src_contig_end = exit_interval.contig_end;
            assoc.aln_src_is_reverse = exit_interval.is_reverse;
            assoc.aln_tgt_contig = entry_interval.contig;
            assoc.aln_tgt_read_start = entry_interval.read_start_est;
            assoc.aln_tgt_read_end = entry_interval.read_end_est;
            assoc.aln_tgt_contig_start = entry_interval.contig_start;
            assoc.aln_tgt_contig_end = entry_interval.contig_end;
            assoc.aln_tgt_is_reverse = entry_interval.is_reverse;
            adjacency_reads.push_back(assoc);
        }
    }
    
    // track cluster associations if clusters are loaded
    if (has_clusters) {
        track_cluster_association(seg_src, seg_tgt, side_src, side_tgt, is_adjacent, read_id);
    }
}

void SegMatrix::track_cluster_association(const string& seg_src, const string& seg_tgt,
                                          const string& side_src, const string& side_tgt,
                                          bool is_adjacent, const string& read_id)
{
    auto seg_src_cluster_it = segment_to_cluster.find(seg_src);
    auto seg_tgt_cluster_it = segment_to_cluster.find(seg_tgt);
    
    if (seg_src_cluster_it == segment_to_cluster.end()) {
        massert(false, "segment %s not found in cluster mapping", seg_src.c_str());
    }
    if (seg_tgt_cluster_it == segment_to_cluster.end()) {
        massert(false, "segment %s not found in cluster mapping", seg_tgt.c_str());
    }
    
    string cluster_src = seg_src_cluster_it->second;
    string cluster_tgt = seg_tgt_cluster_it->second;
    
    // determine cluster sides based on segment sides and strand
    // if segment strand is "-", then left of segment is right of cluster
    string cluster_side_src = side_src;
    string cluster_side_tgt = side_tgt;
    
    auto seg_src_strand_it = segment_to_strand.find(seg_src);
    auto seg_tgt_strand_it = segment_to_strand.find(seg_tgt);
    
    if (seg_src_strand_it != segment_to_strand.end() && seg_src_strand_it->second == "-") {
        cluster_side_src = (side_src == "left") ? "right" : "left";
    }
    if (seg_tgt_strand_it != segment_to_strand.end() && seg_tgt_strand_it->second == "-") {
        cluster_side_tgt = (side_tgt == "left") ? "right" : "left";
    }
    
    auto cluster_key = canonize_key(cluster_src, cluster_tgt, cluster_side_src, cluster_side_tgt);
    
    // track unique reads for cluster associations
    cluster_reach_reads[cluster_key].insert(read_id);
    if (is_adjacent) {
        cluster_adjacency_reads[cluster_key].insert(read_id);
    }
    
    // track associated clusters in cluster stats
    ClusterStats& cluster_src_stats = cluster_stats[cluster_src];
    ClusterStats& cluster_tgt_stats = cluster_stats[cluster_tgt];
    
    if (cluster_side_src == "left") {
        cluster_src_stats.associated_clusters_left.insert(cluster_tgt);
    } else {
        cluster_src_stats.associated_clusters_right.insert(cluster_tgt);
    }
    if (cluster_side_tgt == "left") {
        cluster_tgt_stats.associated_clusters_left.insert(cluster_src);
    } else {
        cluster_tgt_stats.associated_clusters_right.insert(cluster_src);
    }
}

bool SegMatrix::process_interval_side(const ReadInterval& interval,
                                     const string& direction,
                                     bool debug_mode,
                                     unordered_map<string, SegmentReadContribution>& read_contrib,
                                     string& side_out,
                                     double& mutation_percent_out)
{
    if (direction == "entry") {
        side_out = interval.is_reverse ? "right" : "left";
        mutation_percent_out = interval.is_reverse ?
            interval.right_mutation_percent : interval.left_mutation_percent;
    } else {
        side_out = interval.is_reverse ? "left" : "right";
        mutation_percent_out = interval.is_reverse ?
            interval.left_mutation_percent : interval.right_mutation_percent;
    }
    
    if (should_skip_mutation_threshold(mutation_percent_out, interval, direction, debug_mode)) {
        return false;
    }
    
    if (debug_mode) {
        if (direction == "entry") {
            print_entry_detected(interval, mutation_percent_out);
        } else {
            print_exit_detected(interval, mutation_percent_out);
        }
    }
    
    SegmentReadContribution& contrib = read_contrib[interval.segment_id];
    if (side_out == "left") {
        contrib.total_left = true;
    } else {
        contrib.total_right = true;
    }
    
    return true;
}

void SegMatrix::process_interval_sequence(const vector<ReadInterval>& intervals, bool debug_mode, const string& read_id)
{
    if (intervals.empty()) {
        return;
    }
    
    if (debug_mode) {
        cout << "processing interval sequence..." << endl;
    }
    
    std::unordered_map<std::string, SegmentReadContribution> read_contrib;
    
    // track all "going out" intervals seen so far
    vector<const ReadInterval*> going_out_intervals;
    
    // single loop through intervals
    for (const ReadInterval& interval : intervals) {
        if (interval.is_short_segment) {
            if (debug_mode) {
                print_short_segment_skipped(interval);
            }
            continue;
        }
 
        // handle "going in" intervals first
        if (is_pointing_in(interval)) {
            string entry_side;
            double entry_mutation_percent;
            if (process_interval_side(interval, "entry", debug_mode, read_contrib, entry_side, entry_mutation_percent)) {
                if (!interval.is_short_segment) {
                    for (const ReadInterval* exit_interval : going_out_intervals) {
                        if (!exit_interval->is_short_segment) {
                            create_association(*exit_interval, interval, debug_mode, read_contrib, read_id);
                        }
                    }
                }
            }
        }
 
        // handle "going out" intervals after processing entries
        if (is_pointing_out(interval) && !interval.is_short_segment) {
            string exit_side;
            double exit_mutation_percent;
            if (process_interval_side(interval, "exit", debug_mode, read_contrib, exit_side, exit_mutation_percent)) {
                going_out_intervals.push_back(&interval);
            }
        }
    }
    
    for (const auto& kv : read_contrib) {
        const string& seg_id = kv.first;
        const SegmentReadContribution& contrib = kv.second;
        SegmentStats& stats = segment_stats[seg_id];
        if (contrib.total_left) {
            stats.total_read_count_left++;
        }
        if (contrib.total_right) {
            stats.total_read_count_right++;
        }
        if (contrib.associated_left) {
            stats.associated_read_count_left++;
        }
        if (contrib.associated_right) {
            stats.associated_read_count_right++;
        }
        
        // track cluster read contributions (unique reads per cluster side)
        if (has_clusters) {
            track_cluster_read_contribution(seg_id, contrib, read_id);
        }
    }
    
    if (debug_mode) {
        cout << "finished processing interval sequence" << endl;
    }
}

void SegMatrix::track_cluster_read_contribution(const string& seg_id, 
                                                 const SegmentReadContribution& contrib,
                                                 const string& read_id)
{
    auto cluster_it = segment_to_cluster.find(seg_id);
    if (cluster_it == segment_to_cluster.end()) {
        return;
    }
    
    const string& cluster_id = cluster_it->second;
    auto strand_it = segment_to_strand.find(seg_id);
    string strand = (strand_it != segment_to_strand.end()) ? strand_it->second : "+";
    
    ClusterStats& stats = cluster_stats[cluster_id];
    
    // determine cluster sides based on segment sides and strand
    if (contrib.total_left) {
        string cluster_side = (strand == "-") ? "right" : "left";
        if (cluster_side == "left") {
            stats.total_reads_left.insert(read_id);
        } else {
            stats.total_reads_right.insert(read_id);
        }
    }
    if (contrib.total_right) {
        string cluster_side = (strand == "-") ? "left" : "right";
        if (cluster_side == "left") {
            stats.total_reads_left.insert(read_id);
        } else {
            stats.total_reads_right.insert(read_id);
        }
    }
    if (contrib.associated_left) {
        string cluster_side = (strand == "-") ? "right" : "left";
        if (cluster_side == "left") {
            stats.associated_reads_left.insert(read_id);
        } else {
            stats.associated_reads_right.insert(read_id);
        }
    }
    if (contrib.associated_right) {
        string cluster_side = (strand == "-") ? "left" : "right";
        if (cluster_side == "left") {
            stats.associated_reads_left.insert(read_id);
        } else {
            stats.associated_reads_right.insert(read_id);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// process_read
////////////////////////////////////////////////////////////////////////////////////////////////////

void SegMatrix::process_read(uint32_t read_index)
{
    // get all alignments for this read, sorted by read_start
    vector<const Alignment*> alignments = get_sorted_alignments_for_read(read_index);
    
    if (alignments.empty()) {
        return;
    }
    
    string read_id = store.get_read_id(read_index);
    
    // save original alignments for debug output
    vector<const Alignment*> alignments_before = alignments;
    
    // filter overlapping alignments
    alignments = get_parsimony_read_alignments(alignments);
    
    // build sequence of intervals
    vector<ReadInterval> all_intervals;
    
    for (const Alignment* aln : alignments) {
        // intersect with segments
        vector<ReadInterval> intervals = intersect_alignment_with_segments(*aln);
        
        // compute side coverage and mutations for each interval
        for (ReadInterval& interval : intervals) {
            compute_side_coverage_and_mutations(interval, *aln);
            all_intervals.push_back(interval);
        }
    }
    
    // determine debug mode
    string debug_type = "read";  // "off", "few", or "read"
    bool debug_mode = false;
    if (debug_type == "read") {
        debug_mode = (read_id == "m84221_250708_233611_s3/221643819/ccs");
    } else if (debug_type == "few") {
        debug_mode = (debug_read_count < 10) && read_covers_seam(alignments);
    }
    
    // debug: print read info, alignments, and intervals
    if (debug_mode) {
        print_read_debug_info(read_index, read_id, alignments_before, alignments, all_intervals);
    }
    
    // process interval sequence
    process_interval_sequence(all_intervals, debug_mode, read_id);
    
    // increment debug counter only for reads with multiple intervals that we debugged
    if (debug_mode) {
        debug_read_count++;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// compute (main entry point)
////////////////////////////////////////////////////////////////////////////////////////////////////

void SegMatrix::compute(const string& ifn_segments,
                       const string& odir,
                       double max_mutation_percent_param,
                       uint32_t max_adjacency_distance_param,
                       uint32_t max_margin_param,
                       int min_indel_length_param,
                       uint32_t side_length_param,
                       uint32_t side_margin_param,
                       bool output_read_details,
                       const string& ifn_segment_clusters)
{
    max_mutation_percent = max_mutation_percent_param;
    max_adjacency_distance = max_adjacency_distance_param;
    max_margin = max_margin_param;
    min_indel_length = min_indel_length_param;
    side_length = side_length_param;
    side_margin = side_margin_param;
    min_segment_length = side_length + side_margin;
    massert(side_length > 0, "side_length must be positive");
    output_read_details_flag = output_read_details;
    
    segment_stats.clear();
    adjacency_matrix.clear();
    reach_matrix.clear();
    adjacency_reads.clear();
    reach_reads.clear();
    cluster_stats.clear();
    cluster_adjacency_reads.clear();
    cluster_reach_reads.clear();
    debug_read_count = 0;
    
    load_segments(ifn_segments);
    load_cluster_mapping(ifn_segment_clusters);
    
    // validate all segments have clusters if clusters are loaded
    if (has_clusters) {
        for (const auto& seg : segments) {
            if (segment_to_cluster.find(seg.id) == segment_to_cluster.end()) {
                massert(false, "segment %s does not have a cluster", seg.id.c_str());
            }
        }
    }
    
    build_contig_to_segments_map();
    load_library();
    
    // initialize local align if needed (hardcoded: clip_mode=all)
    init_local_align_if_needed(store, ClipMode::ALL);
    
    // process all reads
    const vector<Read>& reads = store.get_reads();
    cout << "processing " << reads.size() << " reads..." << endl;
    
    uint32_t processed = 0;
    for (uint32_t i = 0; i < reads.size(); ++i) {
        process_read(i);
        processed++;
        
        if (processed % 100000 == 0) {
            cout << "processed " << processed << "/" << reads.size() << " reads" << endl;
        }
    }
    
    cout << "finished processing reads" << endl;
    cout << "adjacency matrix entries: " << adjacency_matrix.size() << endl;
    cout << "reach matrix entries: " << reach_matrix.size() << endl;
    
    write_matrices(odir);
    write_segment_summary(odir);
    
    if (output_read_details_flag) {
        write_read_details(odir);
    }
    
    if (has_clusters) {
        compute_cluster_matrices();
        write_cluster_matrices(odir);
        write_cluster_summary(odir);
    }
}

bool SegMatrix::read_covers_seam(const vector<const Alignment*>& alignments) const
{
    if (alignments.size() != 2) {
        return false;
    }
 
    size_t contig_index = alignments[0]->contig_index;
    if (alignments[1]->contig_index != contig_index) {
        return false;
    }
 
    const auto& contigs = store.get_contigs();
    massert(contig_index < contigs.size(), "invalid contig_index %zu", contig_index);
    uint32_t contig_length = contigs[contig_index].length;
    if (contig_length <= 100000) {
        return false;
    }
 
    bool has_start_cover = false;
    bool has_end_cover = false;

    for (const Alignment* aln : alignments) {
        if (!has_start_cover && aln->contig_start <= 10) {
            has_start_cover = true;
        }
        if (!has_end_cover && aln->contig_end >= (contig_length - 10)) {
            has_end_cover = true;
        }
    }

    return has_start_cover && has_end_cover;
}

std::tuple<std::string, std::string, std::string, std::string> 
SegMatrix::canonize_key(const std::string& seg_src, const std::string& seg_tgt,
                       const std::string& side_src, const std::string& side_tgt) const
{
    // compare pairs lexicographically
    if (std::make_pair(seg_src, side_src) > std::make_pair(seg_tgt, side_tgt)) {
        return std::make_tuple(seg_tgt, seg_src, side_tgt, side_src);
    }
    return std::make_tuple(seg_src, seg_tgt, side_src, side_tgt);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// cluster matrix computation
////////////////////////////////////////////////////////////////////////////////////////////////////

void SegMatrix::compute_cluster_matrices()
{
    // cluster matrices are already computed during read processing
    // (cluster_adjacency_reads and cluster_reach_reads track unique reads)
    // this function is a placeholder for any additional computation needed
    cout << "cluster adjacency pairs: " << cluster_adjacency_reads.size() << endl;
    cout << "cluster reach pairs: " << cluster_reach_reads.size() << endl;
}

void SegMatrix::write_cluster_matrices(const string& odir)
{
    // write cluster adjacency matrix
    string adjacency_file = odir + "/cluster_adjacency_matrix.txt";
    ofstream out_adj(adjacency_file);
    massert(out_adj.is_open(), "could not open file %s", adjacency_file.c_str());
    
    out_adj << "cluster_src\tcluster_tgt\tside_src\tside_tgt\ttotal_read_count\tassociated_read_count\tcount" << endl;
    
    for (const auto& entry : cluster_adjacency_reads) {
        const auto& key = entry.first;
        const unordered_set<string>& reads = entry.second;
        uint32_t count = reads.size();
        
        const string& cluster_src = get<0>(key);
        const string& cluster_tgt = get<1>(key);
        const string& side_src = get<2>(key);
        const string& side_tgt = get<3>(key);
        
        uint64_t total = 0;
        uint64_t assoc = 0;
        auto stats_it = cluster_stats.find(cluster_src);
        if (stats_it != cluster_stats.end()) {
            const ClusterStats& stats = stats_it->second;
            total = (side_src == "left") ? stats.total_reads_left.size() : stats.total_reads_right.size();
            assoc = (side_src == "left") ? stats.associated_reads_left.size() : stats.associated_reads_right.size();
        }
        
        out_adj << cluster_src << "\t" << cluster_tgt << "\t"
                << side_src << "\t" << side_tgt << "\t"
                << total << "\t" << assoc << "\t"
                << count << endl;
    }
    
    // write reverse entries (for symmetry)
    for (const auto& entry : cluster_adjacency_reads) {
        const auto& key = entry.first;
        const unordered_set<string>& reads = entry.second;
        uint32_t count = reads.size();
        
        const string& cluster_src = get<0>(key);
        const string& cluster_tgt = get<1>(key);
        if (cluster_src == cluster_tgt) {
            continue;
        }
        const string& side_src = get<2>(key);
        const string& side_tgt = get<3>(key);
        
        uint64_t total = 0;
        uint64_t assoc = 0;
        auto stats_it = cluster_stats.find(cluster_tgt);
        if (stats_it != cluster_stats.end()) {
            const ClusterStats& stats = stats_it->second;
            total = (side_tgt == "left") ? stats.total_reads_left.size() : stats.total_reads_right.size();
            assoc = (side_tgt == "left") ? stats.associated_reads_left.size() : stats.associated_reads_right.size();
        }
        
        out_adj << cluster_tgt << "\t" << cluster_src << "\t"
                << side_tgt << "\t" << side_src << "\t"
                << total << "\t" << assoc << "\t"
                << count << endl;
    }
    
    out_adj.close();
    cout << "wrote cluster adjacency matrix: " << adjacency_file << " (" 
         << cluster_adjacency_reads.size() << " entries)" << endl;
    
    // write cluster reach matrix
    string reach_file = odir + "/cluster_reach_matrix.txt";
    ofstream out_reach(reach_file);
    massert(out_reach.is_open(), "could not open file %s", reach_file.c_str());
    
    out_reach << "cluster_src\tcluster_tgt\tside_src\tside_tgt\ttotal_read_count\tassociated_read_count\tcount" << endl;
    
    for (const auto& entry : cluster_reach_reads) {
        const auto& key = entry.first;
        const unordered_set<string>& reads = entry.second;
        uint32_t count = reads.size();
        
        const string& cluster_src = get<0>(key);
        const string& cluster_tgt = get<1>(key);
        const string& side_src = get<2>(key);
        const string& side_tgt = get<3>(key);
        
        uint64_t total = 0;
        uint64_t assoc = 0;
        auto stats_it = cluster_stats.find(cluster_src);
        if (stats_it != cluster_stats.end()) {
            const ClusterStats& stats = stats_it->second;
            total = (side_src == "left") ? stats.total_reads_left.size() : stats.total_reads_right.size();
            assoc = (side_src == "left") ? stats.associated_reads_left.size() : stats.associated_reads_right.size();
        }
        
        out_reach << cluster_src << "\t" << cluster_tgt << "\t"
                  << side_src << "\t" << side_tgt << "\t"
                  << total << "\t" << assoc << "\t"
                  << count << endl;
    }
    
    // write reverse entries (for symmetry)
    for (const auto& entry : cluster_reach_reads) {
        const auto& key = entry.first;
        const unordered_set<string>& reads = entry.second;
        uint32_t count = reads.size();
        
        const string& cluster_src = get<0>(key);
        const string& cluster_tgt = get<1>(key);
        if (cluster_src == cluster_tgt) {
            continue;
        }
        const string& side_src = get<2>(key);
        const string& side_tgt = get<3>(key);
        
        uint64_t total = 0;
        uint64_t assoc = 0;
        auto stats_it = cluster_stats.find(cluster_tgt);
        if (stats_it != cluster_stats.end()) {
            const ClusterStats& stats = stats_it->second;
            total = (side_tgt == "left") ? stats.total_reads_left.size() : stats.total_reads_right.size();
            assoc = (side_tgt == "left") ? stats.associated_reads_left.size() : stats.associated_reads_right.size();
        }
        
        out_reach << cluster_tgt << "\t" << cluster_src << "\t"
                  << side_tgt << "\t" << side_src << "\t"
                  << total << "\t" << assoc << "\t"
                  << count << endl;
    }
    
    out_reach.close();
    cout << "wrote cluster reach matrix: " << reach_file << " (" 
         << cluster_reach_reads.size() << " entries)" << endl;
}

void SegMatrix::write_cluster_summary(const string& odir)
{
    string summary_file = odir + "/cluster_summary.txt";
    
    ofstream out(summary_file);
    massert(out.is_open(), "could not open file %s", summary_file.c_str());
    
    out << "cluster\ttotal_read_count_left\ttotal_read_count_right\tassociated_read_count_left\tassociated_read_count_right\tassociated_cluster_count_left\tassociated_cluster_count_right" << endl;
    
    // collect all unique clusters from segment_to_cluster to ensure we output all of them
    unordered_set<string> all_clusters;
    for (const auto& kv : segment_to_cluster) {
        all_clusters.insert(kv.second);
    }
    
    // iterate through all clusters to ensure we output all of them
    for (const string& cluster_id : all_clusters) {
        auto stats_it = cluster_stats.find(cluster_id);
        
        uint64_t total_left = 0;
        uint64_t total_right = 0;
        uint64_t assoc_left = 0;
        uint64_t assoc_right = 0;
        uint64_t assoc_cluster_left = 0;
        uint64_t assoc_cluster_right = 0;
        
        if (stats_it != cluster_stats.end()) {
            const ClusterStats& stats = stats_it->second;
            total_left = stats.total_reads_left.size();
            total_right = stats.total_reads_right.size();
            assoc_left = stats.associated_reads_left.size();
            assoc_right = stats.associated_reads_right.size();
            assoc_cluster_left = stats.associated_clusters_left.size();
            assoc_cluster_right = stats.associated_clusters_right.size();
        }
        
        out << cluster_id << "\t" << total_left << "\t" << total_right
            << "\t" << assoc_left << "\t" << assoc_right
            << "\t" << assoc_cluster_left << "\t" << assoc_cluster_right
            << endl;
    }
    
    out.close();
    cout << "wrote cluster summary: " << summary_file << " (" 
         << all_clusters.size() << " clusters)" << endl;
}

