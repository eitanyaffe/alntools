#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Params.h"
#include "alignment_store.h"
#include "utils.h"

using namespace std;

struct SegmentRecord {
    string raw_line;         // original line verbatim (for output)
    string contig;
    uint32_t start;          // 1-based
    uint32_t end;            // 1-based inclusive
    uint32_t length;
    int contig_idx;          // index into store (-1 if not found)
};

// read segment table header and rows, record column indices for contig/start/end
static void read_segment_table(const string& filename,
    vector<string>& header_fields,
    vector<SegmentRecord>& records)
{
    ifstream in(filename);
    massert(in.is_open(), "failed to open segment table: %s", filename.c_str());

    string line;
    massert(bool(getline(in, line)), "segment table is empty: %s", filename.c_str());

    // parse header
    {
        istringstream ss(line);
        string field;
        while (getline(ss, field, '\t'))
            header_fields.push_back(field);
    }

    // find column indices
    int col_contig = -1, col_start = -1, col_end = -1;
    for (int i = 0; i < (int)header_fields.size(); ++i) {
        if (header_fields[i] == "contig") col_contig = i;
        else if (header_fields[i] == "start") col_start = i;
        else if (header_fields[i] == "end")   col_end = i;
    }
    massert(col_contig >= 0, "segment table missing 'contig' column: %s", filename.c_str());
    massert(col_start  >= 0, "segment table missing 'start' column: %s", filename.c_str());
    massert(col_end    >= 0, "segment table missing 'end' column: %s", filename.c_str());

    while (getline(in, line)) {
        if (line.empty()) continue;
        SegmentRecord rec;
        rec.raw_line = line;

        // parse only what we need (contig, start, end) by splitting on tab
        vector<string> fields;
        {
            istringstream ss(line);
            string field;
            while (getline(ss, field, '\t'))
                fields.push_back(field);
        }
        massert((int)fields.size() > col_end,
            "segment table row has too few columns in: %s", filename.c_str());
        rec.contig = fields[col_contig];
        rec.start  = (uint32_t)stoul(fields[col_start]);
        rec.end    = (uint32_t)stoul(fields[col_end]);
        rec.length = rec.end - rec.start + 1;
        rec.contig_idx = -1;
        records.push_back(move(rec));
    }
}

// compute x_coverage for each segment: fraction of bases covered by >=1 passing alignment
static vector<double> compute_coverage(const AlignmentStore& store,
    vector<SegmentRecord>& records,
    ClipMode clip_mode,
    int clip_margin,
    double min_mutations_percent,
    double max_mutations_percent,
    int min_alignment_length,
    int max_alignment_length)
{
    // resolve contig indices once
    for (auto& rec : records) {
        if (store.has_contig_id(rec.contig))
            rec.contig_idx = (int)store.get_contig_index(rec.contig);
    }

    vector<double> coverages(records.size(), 0.0);

    for (size_t i = 0; i < records.size(); ++i) {
        const SegmentRecord& rec = records[i];
        if (rec.contig_idx < 0 || rec.length == 0) {
            coverages[i] = 0.0;
            continue;
        }

        // segment in 0-based half-open coordinates
        uint32_t seg_start = rec.start - 1;
        uint32_t seg_end   = rec.end;       // exclusive

        // bitmap over segment positions
        vector<bool> bitmap(rec.length, false);

        Interval interval(rec.contig, seg_start, seg_end);
        auto alns = store.get_alignments_intersecting_interval(interval);

        for (const Alignment& aln : alns) {
            if (!passes_alignment_filter(aln, store, clip_mode, clip_margin,
                    min_mutations_percent, max_mutations_percent,
                    min_alignment_length, max_alignment_length))
                continue;

            // clamp alignment to segment
            uint32_t ov_start = max(aln.contig_start, seg_start);
            uint32_t ov_end   = min(aln.contig_end,   seg_end);
            if (ov_start >= ov_end) continue;

            // mark covered positions relative to segment start
            for (uint32_t pos = ov_start; pos < ov_end; ++pos)
                bitmap[pos - seg_start] = true;
        }

        uint32_t covered = 0;
        for (bool b : bitmap)
            if (b) ++covered;
        coverages[i] = (double)covered / rec.length;
    }

    return coverages;
}

void cov_intervals_command(const string& aln_file,
    const string& ifn_segments,
    const string& ofn,
    ClipMode clip_mode,
    int clip_margin,
    double min_mutations_percent,
    double max_mutations_percent,
    int min_alignment_length,
    int max_alignment_length,
    int min_indel_length)
{
    cout << "loading alignment file " << aln_file << endl;
    AlignmentStore store;
    store.load(aln_file);
    store.count_short_indels(min_indel_length);

    cout << "reading segment table " << ifn_segments << endl;
    vector<string> header_fields;
    vector<SegmentRecord> records;
    read_segment_table(ifn_segments, header_fields, records);
    cout << "segments: " << records.size() << endl;

    vector<double> coverages = compute_coverage(store, records,
        clip_mode, clip_margin,
        min_mutations_percent, max_mutations_percent,
        min_alignment_length, max_alignment_length);

    cout << "writing output to " << ofn << endl;
    ofstream out(ofn);
    massert(out.is_open(), "failed to open output file: %s", ofn.c_str());

    // header with appended x_coverage column
    for (size_t i = 0; i < header_fields.size(); ++i) {
        if (i > 0) out << "\t";
        out << header_fields[i];
    }
    out << "\tx_coverage\n";

    for (size_t i = 0; i < records.size(); ++i) {
        out << records[i].raw_line << "\t" << coverages[i] << "\n";
    }

    cout << "cov_intervals done" << endl;
}

void cov_intervals_params(const char* name, int argc, char** argv, Parameters& params)
{
    params.add_parser("ifn_aln",               new ParserFilename("input ALN file"),                              true);
    params.add_parser("ifn_segments",          new ParserFilename("input segment table (contig, start, end, ...)"), true);
    params.add_parser("ofn",                   new ParserFilename("output file"),                                  true);
    params.add_parser("clip_mode",             new ParserString("clip mode", "complete"),                          false);
    params.add_parser("clip_margin",           new ParserInteger("clip margin in bases", 10),                      false);
    params.add_parser("min_mutations_percent", new ParserDouble("min mutation percent", 0.0),                      false);
    params.add_parser("max_mutations_percent", new ParserDouble("max mutation percent", 0.1),                      false);
    params.add_parser("min_alignment_length",  new ParserInteger("min alignment length (0=none)", 1000),           false);
    params.add_parser("max_alignment_length",  new ParserInteger("max alignment length (0=none)", 0),              false);
    params.add_parser("min_indel_length",      new ParserInteger("min indel length", 3),                           false);

    if (argc == 1) {
        params.usage(name);
        exit(1);
    }
    params.read(argc, argv);
    params.parse();
    params.verify_mandatory();
    params.print(cout);
}

int cov_intervals_main(const char* name, int argc, char** argv)
{
    Parameters params;
    cov_intervals_params(name, argc, argv, params);

    string ifn_aln      = params.get_string("ifn_aln");
    string ifn_segments = params.get_string("ifn_segments");
    string ofn          = params.get_string("ofn");

    ClipMode clip_mode  = string_to_clip_mode(params.get_string("clip_mode"));
    int    clip_margin  = params.get_int("clip_margin");
    double min_mut      = params.get_double("min_mutations_percent");
    double max_mut      = params.get_double("max_mutations_percent");
    int    min_aln_len  = params.get_int("min_alignment_length");
    int    max_aln_len  = params.get_int("max_alignment_length");
    int    min_indel    = params.get_int("min_indel_length");

    cov_intervals_command(ifn_aln, ifn_segments, ofn,
        clip_mode, clip_margin, min_mut, max_mut,
        min_aln_len, max_aln_len, min_indel);

    return 0;
}
