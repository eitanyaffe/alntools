#include "GetReadIds.h"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

GetReadIds::GetReadIds() {}

void GetReadIds::load_segments(const string& ifn_segments)
{
    cout << "reading segment table: " << ifn_segments << endl;
    ifstream file(ifn_segments);
    massert(file.is_open(), "could not open file %s", ifn_segments.c_str());

    string line;
    int contig_ind = -1, start_ind = -1, end_ind = -1, bin_id_ind = -1;

    if (getline(file, line)) {
        istringstream header_iss(line);
        string header;
        int col_index = 0;
        while (getline(header_iss, header, '\t')) {
            if (header == "contig") contig_ind = col_index;
            else if (header == "start") start_ind = col_index;
            else if (header == "end") end_ind = col_index;
            else if (header == "bin_id") bin_id_ind = col_index;
            col_index++;
        }
        massert(contig_ind >= 0, "column 'contig' not found in %s", ifn_segments.c_str());
        massert(start_ind >= 0, "column 'start' not found in %s", ifn_segments.c_str());
        massert(end_ind >= 0, "column 'end' not found in %s", ifn_segments.c_str());
        massert(bin_id_ind >= 0, "column 'bin_id' not found in %s", ifn_segments.c_str());
    }

    while (getline(file, line)) {
        if (line.empty()) continue;

        istringstream iss(line);
        vector<string> fields;
        string field;
        while (getline(iss, field, '\t'))
            fields.push_back(field);

        if (fields.empty()) continue;

        ReadSegment seg;
        seg.contig = fields[contig_ind];
        seg.start  = atoi(fields[start_ind].c_str());
        seg.end    = atoi(fields[end_ind].c_str());
        seg.bin_id = fields[bin_id_ind];
        segments.push_back(seg);
    }

    cout << "loaded " << segments.size() << " segments" << endl;
}

void GetReadIds::collect_read_ids(AlignmentStore& store,
                                   unordered_map<uint32_t, ReadBinEntry>& read_map) const
{
    store.count_short_indels(min_indel_length);
    init_local_align_if_needed(store, clip_mode);

    for (const auto& seg : segments) {
        // convert 1-based inclusive to 0-based half-open
        uint32_t start_0based = seg.start - 1;
        uint32_t end_0based   = seg.end;

        Interval interval(seg.contig, start_0based, end_0based);
        auto alignments = store.get_alignments_intersecting_interval(interval);

        for (const auto& aln_ref : alignments) {
            const Alignment& aln = aln_ref.get();

            if (!passes_alignment_filter(aln, store, clip_mode, clip_margin,
                                         min_mutations_percent, max_mutations_percent,
                                         min_alignment_length, max_alignment_length)) {
                continue;
            }

            // intersection length in contig coordinates
            uint32_t isect_start = max(start_0based, aln.contig_start);
            uint32_t isect_end   = min(end_0based,   aln.contig_end);
            if (isect_end <= isect_start) continue;
            uint32_t isect_len = isect_end - isect_start;

            uint32_t read_index = aln.read_index;
            auto it = read_map.find(read_index);
            if (it == read_map.end() || isect_len > it->second.length) {
                read_map[read_index] = ReadBinEntry(seg.bin_id, isect_len);
            }
        }
    }
}

void GetReadIds::write_output(const string& ofn,
                               const AlignmentStore& store,
                               const unordered_map<uint32_t, ReadBinEntry>& read_map) const
{
    cout << "writing output: " << ofn << endl;
    ofstream out(ofn.c_str());
    massert(out.is_open(), "could not open file %s", ofn.c_str());

    out << "read_id\tbin_id\n";
    const auto& reads = store.get_reads();
    for (uint32_t i = 0; i < static_cast<uint32_t>(reads.size()); ++i) {
        auto it = read_map.find(i);
        const string& bin_id = (it != read_map.end()) ? it->second.bin_id : "none";
        out << reads[i].id << "\t" << bin_id << "\n";
    }

    out.close();
}

void GetReadIds::compute(const string& ifn_aln,
                          const string& ifn_segments,
                          const string& ofn,
                          ClipMode clip_mode_param,
                          int clip_margin_param,
                          double min_mutations_percent_param,
                          double max_mutations_percent_param,
                          int min_alignment_length_param,
                          int max_alignment_length_param,
                          int min_indel_length_param)
{
    clip_mode             = clip_mode_param;
    clip_margin           = clip_margin_param;
    min_mutations_percent = min_mutations_percent_param;
    max_mutations_percent = max_mutations_percent_param;
    min_alignment_length  = min_alignment_length_param;
    max_alignment_length  = max_alignment_length_param;
    min_indel_length      = min_indel_length_param;

    load_segments(ifn_segments);

    cout << "loading alignment file: " << ifn_aln << endl;
    AlignmentStore store;
    store.load(ifn_aln);

    unordered_map<uint32_t, ReadBinEntry> read_map;
    collect_read_ids(store, read_map);
    cout << "assigned bin to " << read_map.size() << " reads" << endl;

    write_output(ofn, store, read_map);
    cout << "get_read_ids complete" << endl;
}
