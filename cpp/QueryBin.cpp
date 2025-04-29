#include "QueryBin.h"
#include <algorithm> // For std::min/max
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

QueryBin::QueryBin(
    const std::vector<Interval>& intervals,
    const AlignmentStore& store,
    int binsize)
    : intervals(intervals)
    , store(store)
    , binsize(binsize)
{
  if (binsize <= 0) {
    cerr << "error: binsize must be positive." << endl;
    exit(1);
  }
}

void QueryBin::aggregate_data()
{
  bin_results.clear();

  for (const auto& interval : intervals) {
    uint32_t contig_index = store.get_contig_index(interval.contig);

    // Calculate relevant bin range based on original interval
    uint32_t adjusted_start = (interval.start / binsize) * binsize;
    // Handle edge case where end is 0 or interval is empty
    if (interval.end == 0 || interval.start >= interval.end)
      continue;
    uint32_t last_bin_start = ((interval.end - 1) / binsize) * binsize;

    // Initialize Relevant Bins in the map
    for (uint32_t b_start = adjusted_start; b_start <= last_bin_start; b_start += binsize) {
      bin_results.try_emplace({ contig_index, b_start }, BinData());
    }

    // Get alignments overlapping the *original* interval
    std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval);

    for (const auto& alignment_ref : alignments) {
      const auto& aln = alignment_ref.get();

      // Iterate through the relevant bins for this interval
      for (uint32_t b_start = adjusted_start; b_start <= last_bin_start; b_start += binsize) {
        uint32_t b_end = b_start + binsize;

        // Calculate overlap considering alignment, bin, AND original interval boundaries
        uint32_t effective_start = std::max({ aln.contig_start, b_start, interval.start });
        uint32_t effective_end = std::min({ aln.contig_end, b_end, interval.end });

        int overlap_length = (effective_end > effective_start) ? (effective_end - effective_start) : 0;

        if (overlap_length > 0) {
          auto it = bin_results.find({ contig_index, b_start });
          // We should always find it because we pre-populated
          if (it != bin_results.end()) {
            // Using int now, check potential overflow? (unlikely for overlap_length)
            it->second.sequenced_basepairs += overlap_length;
          } else {
            // This case indicates a logic error in initialization or calculation
            cerr << "error: bin " << b_start << " on contig " << interval.contig
                 << " should have been initialized but wasn't." << endl;
          }
        }
      }

      // Process mutations
      for (const auto& mutation : aln.mutations) {
        uint32_t mutation_contig_pos = aln.contig_start + mutation.position;

        // Ignore mutations outside the original interval
        if (mutation_contig_pos < interval.start || mutation_contig_pos >= interval.end) {
          continue;
        }

        uint32_t mutation_bin_start = (mutation_contig_pos / binsize) * binsize;

        // Find the bin in our map (it must be relevant if pos is within interval)
        auto it = bin_results.find({ contig_index, mutation_bin_start });
        if (it != bin_results.end()) {
          it->second.mutation_count++;
        } else {
          // Logic error if mutation is inside interval but bin wasn't initialized.
          // Could happen if interval.end=0 edge case calculation was wrong.
          cerr << "error: bin " << mutation_bin_start << " on contig " << interval.contig
               << " should have been initialized but wasn't." << endl;
        }
      }
    }
  }
}

void QueryBin::generate_output_rows()
{
  output_rows.clear();

  for (const auto& entry : bin_results) {
    const auto& key = entry.first;
    const auto& data = entry.second;

    uint32_t contig_index = key.first;
    uint32_t bin_start = key.second;
    uint32_t bin_end = bin_start + binsize; // Bin end is standard
    int bin_length = binsize; // Bin length is standard

    string contig_id = store.get_contig_id(contig_index);

    output_rows.push_back({ contig_id, bin_start, bin_end, bin_length,
        data.sequenced_basepairs, data.mutation_count });
  }
}

// New function to write the generated rows to a file
void QueryBin::write_to_csv(const std::string& ofn_prefix)
{
  string filename = ofn_prefix + "_bins.tsv";
  cout << "writing bin data rows to " << filename << endl;
  ofstream ofs(filename);

  if (!ofs.is_open()) {
    // Consider throwing an exception instead of exiting
    cerr << "error: could not open file " << filename << endl;
    exit(1);
  }

  // Write header - removed coverage column previously present
  ofs << "contig\tbin_start\tbin_end\tbin_length\tsequenced_bp\tmutation_count\n";

  for (const auto& row : output_rows) {
    ofs << row.contig << "\t"
        << row.bin_start << "\t"
        << row.bin_end << "\t"
        << row.bin_length << "\t"
        << row.sequenced_basepairs << "\t"
        << row.mutation_count << "\n";
  }

  ofs.close();
}

void QueryBin::execute()
{
  aggregate_data();
  generate_output_rows();
}
