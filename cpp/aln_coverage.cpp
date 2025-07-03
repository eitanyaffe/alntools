#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include "Params.h"
#include "alignment_store.h"
#include "utils.h"

using namespace std;

struct IntervalInfo {
  string id;
  uint32_t start; // 1-based
  uint32_t end; // 1-based (inclusive)
  uint32_t length;
};

struct ReadCoverageStats {
  unordered_set<uint32_t> aligned_reads;
  uint64_t non_aligned_read_bp = 0;
  vector<IntervalInfo> non_aligned_read_intervals;
};

struct ContigCoverageStats {
  uint64_t non_aligned_contig_bp = 0;
  vector<IntervalInfo> non_aligned_contig_intervals;
};

// extract unaligned intervals from coverage bitmap
vector<IntervalInfo> extract_unaligned_intervals(const string& id, const vector<bool>& coverage)
{
  vector<IntervalInfo> intervals;
  bool in_unaligned_region = false;
  uint32_t start = 0;

  for (uint32_t pos = 0; pos < coverage.size(); ++pos) {
    if (!coverage[pos]) { // unaligned position
      if (!in_unaligned_region) {
        start = pos;
        in_unaligned_region = true;
      }
    } else { // aligned position
      if (in_unaligned_region) {
        // convert to 1-based coordinates for output
        uint32_t start_1based = start + 1;
        uint32_t end_1based = pos; // pos is exclusive, so end is pos
        intervals.push_back({ id, start_1based, end_1based, pos - start });
        in_unaligned_region = false;
      }
    }
  }

  // handle case where sequence ends with unaligned region
  if (in_unaligned_region) {
    uint32_t start_1based = start + 1;
    uint32_t end_1based = (uint32_t)coverage.size();
    intervals.push_back({ id, start_1based, end_1based, (uint32_t)coverage.size() - start });
  }

  return intervals;
}

// calculate N50 of contig lengths
uint64_t calculate_contig_n50(const AlignmentStore& store)
{
  vector<uint64_t> lengths;
  for (const auto& contig : store.get_contigs()) {
    lengths.push_back(contig.length);
  }

  // sort in descending order
  sort(lengths.begin(), lengths.end(), greater<uint64_t>());

  // calculate total length
  uint64_t total_length = 0;
  for (uint64_t length : lengths) {
    total_length += length;
  }

  // find N50
  uint64_t cumulative_length = 0;
  uint64_t target_length = total_length / 2;

  for (uint64_t length : lengths) {
    cumulative_length += length;
    if (cumulative_length >= target_length) {
      return length;
    }
  }

  return 0; // should not happen if lengths is not empty
}

// calculate read coverage statistics
ReadCoverageStats calculate_read_coverage_stats(const vector<vector<bool>>& read_coverage,
    const AlignmentStore& store,
    size_t read_count)
{
  ReadCoverageStats stats;

  for (size_t read_idx = 0; read_idx < read_count; ++read_idx) {
    bool has_alignment = false;
    uint32_t unaligned_bp = 0;

    for (bool covered : read_coverage[read_idx]) {
      if (covered) {
        has_alignment = true;
      } else {
        unaligned_bp++;
      }
    }

    if (has_alignment) {
      stats.aligned_reads.insert(read_idx);
    }

    stats.non_aligned_read_bp += unaligned_bp;

    // extract unaligned intervals for this read
    auto intervals = extract_unaligned_intervals(store.get_read_id(read_idx), read_coverage[read_idx]);
    stats.non_aligned_read_intervals.insert(stats.non_aligned_read_intervals.end(), intervals.begin(), intervals.end());
  }

  return stats;
}

// calculate contig coverage statistics
ContigCoverageStats calculate_contig_coverage_stats(const vector<vector<bool>>& contig_coverage,
    const AlignmentStore& store,
    size_t contig_count)
{
  ContigCoverageStats stats;

  for (size_t contig_idx = 0; contig_idx < contig_count; ++contig_idx) {
    uint32_t unaligned_bp = 0;

    for (bool covered : contig_coverage[contig_idx]) {
      if (!covered) {
        unaligned_bp++;
      }
    }

    stats.non_aligned_contig_bp += unaligned_bp;

    // extract unaligned intervals for this contig
    auto intervals = extract_unaligned_intervals(store.get_contig_id(contig_idx), contig_coverage[contig_idx]);
    stats.non_aligned_contig_intervals.insert(stats.non_aligned_contig_intervals.end(), intervals.begin(), intervals.end());
  }

  return stats;
}

// write main coverage statistics table
void write_coverage_stats(const string& output_prefix,
    size_t contig_count, size_t read_count, size_t alignment_count,
    uint64_t total_read_bp, uint64_t total_assembly_bp, uint64_t contig_n50, size_t aligned_count,
    const ReadCoverageStats& read_stats, const ContigCoverageStats& contig_stats)
{
  string coverage_file = output_prefix + "_coverage.txt";
  cout << "writing coverage to " << coverage_file << endl;
  ofstream coverage_out(coverage_file);
  massert(coverage_out.is_open(), "failed to open file for writing: %s", coverage_file.c_str());

  coverage_out << "stat\tvalue\n";
  coverage_out << "contig_count\t" << contig_count << "\n";
  coverage_out << "read_count\t" << read_count << "\n";
  coverage_out << "alignment_count\t" << alignment_count << "\n";
  coverage_out << "total_read_bp\t" << total_read_bp << "\n";
  coverage_out << "total_assembly_bp\t" << total_assembly_bp << "\n";
  coverage_out << "contig_n50\t" << contig_n50 << "\n";
  coverage_out << "aligned_count\t" << aligned_count << "\n";
  coverage_out << "non_aligned_read_bp\t" << read_stats.non_aligned_read_bp << "\n";
  coverage_out << "non_aligned_contig_bp\t" << contig_stats.non_aligned_contig_bp << "\n";

  coverage_out.close();
}

// write non-aligned read intervals
void write_read_intervals(const string& output_prefix, const ReadCoverageStats& read_stats)
{
  string read_intervals_file = output_prefix + "_non_aligned_reads.txt";
  cout << "writing non-aligned read intervals to " << read_intervals_file << endl;
  ofstream read_intervals_out(read_intervals_file);
  massert(read_intervals_out.is_open(), "failed to open file for writing: %s", read_intervals_file.c_str());

  read_intervals_out << "read_id\tstart\tend\tlength\n";
  for (const auto& interval : read_stats.non_aligned_read_intervals) {
    read_intervals_out << interval.id << "\t" << interval.start << "\t"
                       << interval.end << "\t" << interval.length << "\n";
  }

  read_intervals_out.close();
}

// write non-aligned contig intervals
void write_contig_intervals(const string& output_prefix, const ContigCoverageStats& contig_stats)
{
  string contig_intervals_file = output_prefix + "_non_aligned_contigs.txt";
  cout << "writing non-aligned contig intervals to " << contig_intervals_file << endl;
  ofstream contig_intervals_out(contig_intervals_file);
  massert(contig_intervals_out.is_open(), "failed to open file for writing: %s", contig_intervals_file.c_str());

  contig_intervals_out << "contig_id\tstart\tend\tlength\n";
  for (const auto& interval : contig_stats.non_aligned_contig_intervals) {
    contig_intervals_out << interval.id << "\t" << interval.start << "\t"
                         << interval.end << "\t" << interval.length << "\n";
  }

  contig_intervals_out.close();
}

void coverage_command(const string& aln_file, const string& output_prefix)
{
  cout << "loading alignment file " << aln_file << endl;

  AlignmentStore store;
  store.load(aln_file);

  // basic stats (straightforward)
  size_t contig_count = store.get_contigs().size();
  size_t read_count = store.get_reads().size();
  size_t alignment_count = store.get_alignment_count();

  // calculate total bp
  uint64_t total_read_bp = 0;
  for (const auto& read : store.get_reads()) {
    total_read_bp += read.length;
  }

  uint64_t total_assembly_bp = 0;
  for (const auto& contig : store.get_contigs()) {
    total_assembly_bp += contig.length;
  }

  // calculate contig N50
  uint64_t contig_n50 = calculate_contig_n50(store);

  // prepare for bitmap analysis
  vector<vector<bool>> read_coverage(read_count);
  vector<vector<bool>> contig_coverage(contig_count);

  // initialize coverage bitmaps
  for (size_t i = 0; i < read_count; ++i) {
    read_coverage[i].resize(store.get_reads()[i].length, false);
  }
  for (size_t i = 0; i < contig_count; ++i) {
    contig_coverage[i].resize(store.get_contigs()[i].length, false);
  }

  // mark aligned regions in bitmaps
  for (const auto& aln : store.get_alignments()) {
    // mark read coverage
    for (uint32_t pos = aln.read_start; pos < aln.read_end; ++pos) {
      read_coverage[aln.read_index][pos] = true;
    }

    // mark contig coverage
    for (uint32_t pos = aln.contig_start; pos < aln.contig_end; ++pos) {
      contig_coverage[aln.contig_index][pos] = true;
    }
  }

  // calculate read coverage statistics
  ReadCoverageStats read_stats = calculate_read_coverage_stats(read_coverage, store, read_count);

  // calculate contig coverage statistics
  ContigCoverageStats contig_stats = calculate_contig_coverage_stats(contig_coverage, store, contig_count);

  size_t aligned_count = read_stats.aligned_reads.size();

  // write output files
  write_coverage_stats(output_prefix, contig_count, read_count, alignment_count,
      total_read_bp, total_assembly_bp, contig_n50, aligned_count, read_stats, contig_stats);
  write_read_intervals(output_prefix, read_stats);
  write_contig_intervals(output_prefix, contig_stats);

  cout << "coverage calculation complete" << endl;
  cout << "wrote " << read_stats.non_aligned_read_intervals.size() << " non-aligned read intervals" << endl;
  cout << "wrote " << contig_stats.non_aligned_contig_intervals.size() << " non-aligned contig intervals" << endl;
}

void coverage_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("input ALN file"), true);
  params.add_parser("ofn_prefix", new ParserFilename("output prefix for tables"), true);

  if (argc == 1) {
    params.usage(name);
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

int coverage_main(const char* name, int argc, char** argv)
{
  Parameters params;
  coverage_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  string ofn_prefix = params.get_string("ofn_prefix");

  coverage_command(ifn, ofn_prefix);

  return 0;
}