#include "QueryPileup.h"
#include <algorithm> // For std::sort, std::max
#include <cassert> // For assertions
#include <fstream>
#include <iostream>
#include <map>
#include <numeric> // For std::accumulate
#include <string>
#include <vector>

using namespace std;

QueryPileup::QueryPileup(
    const std::vector<Interval>& intervals,
    const std::string& ofn_prefix,
    const AlignmentStore& store,
    PileupReportMode report_mode)
    : intervals(intervals)
    , ofn_prefix(ofn_prefix)
    , store(store)
    , report_mode(report_mode)
{
  // Constructor implementation (basic initialization done via initializer list)
}

// Private helper function to aggregate data into pileup_results
void QueryPileup::aggregate_data()
{
  pileup_results.clear(); // Ensure map is empty before starting

  cout << "aggregating alignment data for pileup..." << endl;

  // Pre-populate pileup_results map for all positions defined by input intervals.
  cout << "initializing positions based on intervals..." << endl;
  for (const auto& interval : intervals) {
    uint32_t contig_index = store.get_contig_index(interval.contig);
    for (uint32_t pos = interval.start; pos < interval.end; ++pos) {
      pileup_results.try_emplace({ contig_index, pos }, PileupData());
    }
  }
  cout << "initialized " << pileup_results.size() << " positions." << endl;

  // Loop through intervals (again) to process alignments associated with them
  cout << "processing alignments..." << endl;
  int processed_alignments = 0;
  for (const auto& interval : intervals) {
    // Get alignments overlapping this interval
    std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval);

    // Loop through alignments
    for (const auto& alignment_ref : alignments) {
      const auto& aln = alignment_ref.get();
      uint32_t contig_index = aln.contig_index;

      // Calculate coverage for relevant positions.
      for (uint32_t pos = aln.contig_start; pos < aln.contig_end; ++pos) {
        auto it = pileup_results.find({ contig_index, pos });
        if (it != pileup_results.end()) {
          it->second.coverage++;
        }
      }

      // Calculate mutation counts for relevant positions.
      for (const auto& mutation : aln.mutations) {
        uint32_t mutation_contig_pos = aln.contig_start + mutation.position;
        auto it = pileup_results.find({ contig_index, mutation_contig_pos });
        if (it != pileup_results.end()) {
          string mut_str = mutation.to_string();
          it->second.mutation_counts[mut_str]++;
        }
      }
      processed_alignments++;
    }
  }
  cout << "processed " << processed_alignments << " alignments." << endl;
}

// Private helper function to populate output_rows from pileup_results
void QueryPileup::generate_output_rows()
{
  output_rows.clear(); // Ensure vector is empty
  cout << "generating output rows..." << endl;

  // Loop through pileup_results map (sorted by contig_idx, position_0based).
  for (const auto& entry : pileup_results) {
    const auto& key = entry.first;
    const auto& data = entry.second;
    uint32_t contig_index = key.first;
    uint32_t pos_0based = key.second;

    // Apply filtering based on report_mode, coverage, and mutation counts.
    if (report_mode == PileupReportMode::COVERED && data.coverage == 0) {
      continue;
    }
    if (report_mode == PileupReportMode::MUTATED && data.mutation_counts.empty()) {
      continue;
    }

    // Calculate position_1based.
    uint32_t position_1based = pos_0based + 1;

    // Get contig_id.
    string contig_id = store.get_contig_id(contig_index);

    // Calculate ref_count.
    int total_mutated_count = 0;
    for (const auto& mc : data.mutation_counts) {
      total_mutated_count += mc.second;
    }
    int ref_count = data.coverage - total_mutated_count;
    assert(ref_count >= 0 && "Reference count cannot be negative");

    // Create and sort vector of observed variants (by count desc, then name asc).
    std::vector<std::pair<std::string, int>> variants(data.mutation_counts.begin(), data.mutation_counts.end());
    std::sort(variants.begin(), variants.end(),
        [](const auto& a, const auto& b) {
          if (a.second != b.second) {
            return a.second > b.second; // Sort by count descending
          }
          return a.first < b.first; // Sort by variant string ascending for ties
        });

    int cumulative_count_for_pos = 0;

    // Loop through sorted variants, calculate cumsum, create output rows.
    for (const auto& variant_pair : variants) {
      const string& mut_str = variant_pair.first;
      int count = variant_pair.second;
      cumulative_count_for_pos += count;
      output_rows.push_back({ contig_id, position_1based, mut_str, count, data.coverage, cumulative_count_for_pos });
    }

    // Create REF row if ref_count > 0 (or if coverage is 0 but mode is ALL).
    if (ref_count > 0 || (data.coverage == 0 && report_mode == PileupReportMode::ALL)) {
      cumulative_count_for_pos += ref_count;
      output_rows.push_back({ contig_id, position_1based, "REF", ref_count, data.coverage, cumulative_count_for_pos });
    }

    assert(cumulative_count_for_pos == data.coverage && "Cumulative count must equal coverage at end of position");
  }
  cout << "generated " << output_rows.size() << " output rows." << endl;
}

// Private helper function to write output_rows to file
void QueryPileup::write_rows_to_file()
{
  string filename = ofn_prefix + "_pileup.tsv";
  cout << "writing pileup data to " << filename << endl;
  ofstream ofs(filename);

  if (!ofs.is_open()) {
    cerr << "error: could not open file " << filename << endl;
    exit(1); // Abort as per instructions
  }

  // Write header
  ofs << "contig\tposition\tvariant\tcount\tcoverage\tcumsum\n";

  for (const auto& row : output_rows) {
    ofs << row.contig << "\t"
        << row.position << "\t"
        << row.variant << "\t"
        << row.count << "\t"
        << row.coverage << "\t"
        << row.cumsum << "\n";
  }

  ofs.close();
  cout << "wrote " << output_rows.size() << " rows to " << filename << endl;
}

// Public method to orchestrate the process
void QueryPileup::write_to_csv()
{
  aggregate_data();
  generate_output_rows(); // Changed from write_data_to_file
  write_rows_to_file(); // New call
}