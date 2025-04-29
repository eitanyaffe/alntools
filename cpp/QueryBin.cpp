#include "QueryBin.h"
#include <algorithm> // For std::min/max
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

QueryBin::QueryBin(
    const std::vector<Interval>& intervals,
    const std::string& ofn_prefix,
    const AlignmentStore& store,
    int binsize,
    bool skip_empty_bins)
    : intervals(intervals)
    , ofn_prefix(ofn_prefix)
    , store(store)
    , binsize(binsize)
    , skip_empty_bins(skip_empty_bins)
{
  if (binsize <= 0) {
    cerr << "error: binsize must be positive." << endl;
    exit(1); // Abort as per instructions
  }
}

// Private helper function to aggregate data into bin_results
void QueryBin::aggregate_data()
{
  bin_results.clear(); // Ensure map is empty before starting

  cout << "aggregating alignment data into bins..." << endl;

  for (const auto& interval : intervals) {
    // Adjust interval boundaries to multiples of binsize
    uint32_t adjusted_start = (interval.start / binsize) * binsize;
    uint32_t adjusted_end = ((interval.end + binsize - 1) / binsize) * binsize;
    if (adjusted_end == interval.end && adjusted_end > adjusted_start) { // Avoid empty range if end is a multiple
      adjusted_end = ((interval.end - 1 + binsize) / binsize) * binsize; // Round end-1 up
    }

    // Ensure adjusted interval is valid
    if (adjusted_start >= adjusted_end) {
      cout << "warning: skipping interval " << interval.contig << ":" << interval.start << "-" << interval.end
           << " because adjusted range [" << adjusted_start << ", " << adjusted_end << ") is empty or invalid." << endl;
      continue;
    }

    Interval adjusted_interval(interval.contig, adjusted_start, adjusted_end);
    uint32_t contig_index = store.get_contig_index(adjusted_interval.contig);

    // Initialize all bins within the adjusted interval
    for (uint32_t current_bin_start = adjusted_start; current_bin_start < adjusted_end; current_bin_start += binsize) {
      // Use try_emplace to avoid overwriting if multiple intervals cover the same bin
      bin_results.try_emplace({ contig_index, current_bin_start }, BinData());
    }

    // Get alignments overlapping the adjusted interval
    std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(adjusted_interval);

    for (const auto& alignment_ref : alignments) {
      const auto& aln = alignment_ref.get();

      // Calculate the range of global bin indices the alignment overlaps
      uint32_t start_bin_index = aln.contig_start / binsize;
      uint32_t end_bin_index = (aln.contig_end - 1) / binsize;

      // Process base pair overlaps
      for (uint32_t bin_index = start_bin_index; bin_index <= end_bin_index; ++bin_index) {
        uint32_t bin_start = bin_index * binsize;
        // Check if this global bin is within our adjusted interval for initialization robustness
        if (bin_start >= adjusted_start && bin_start < adjusted_end) {
          uint32_t bin_end = bin_start + binsize;
          uint32_t overlap_start = std::max(aln.contig_start, bin_start);
          uint32_t overlap_end = std::min(aln.contig_end, bin_end);
          uint32_t overlap_length = (overlap_end > overlap_start) ? (overlap_end - overlap_start) : 0;

          if (overlap_length > 0) {
            // Use find + update to ensure the key exists (initialized above)
            auto it = bin_results.find({ contig_index, bin_start });
            if (it != bin_results.end()) {
              it->second.sequenced_basepairs += overlap_length;
            } // else: Should not happen due to initialization, but safe check
          }
        }
      }

      // Process mutations
      for (const auto& mutation : aln.mutations) {
        uint32_t mutation_contig_pos = aln.contig_start + mutation.position;
        uint32_t mutation_bin_index = mutation_contig_pos / binsize;
        uint32_t mutation_bin_start = mutation_bin_index * binsize;

        // Check if this global bin is within our adjusted interval
        if (mutation_bin_start >= adjusted_start && mutation_bin_start < adjusted_end) {
          // Use find + update
          auto it = bin_results.find({ contig_index, mutation_bin_start });
          if (it != bin_results.end()) {
            it->second.mutation_count++;
          } // else: Should not happen
        }
      }
    }
  }
}

// Private helper function to write aggregated data to the TSV file
void QueryBin::write_data_to_file()
{
  string filename = ofn_prefix + "_bins.tsv";
  cout << "writing bin data (" << bin_results.size() << " bins) to " << filename << endl;
  ofstream ofs(filename);

  if (!ofs.is_open()) {
    cerr << "error: could not open file " << filename << endl;
    exit(1); // Abort as per instructions
  }

  // Write header
  ofs << "contig\tbin_start\tbin_end\tbin_length\tsequenced_bp\tmutation_count\n";

  for (const auto& entry : bin_results) {
    const auto& key = entry.first;
    const auto& data = entry.second;

    uint32_t contig_index = key.first;
    uint32_t bin_start = key.second;
    uint32_t bin_end = bin_start + binsize;

    // Get contig name from store using index
    string contig_id = store.get_contig_id(contig_index);

    // Skip empty bins if skip_empty_bins is true
    if (skip_empty_bins && data.sequenced_basepairs == 0) {
      continue;
    }

    ofs << contig_id << "\t"
        << bin_start << "\t"
        << bin_end << "\t"
        << binsize << "\t"
        << data.sequenced_basepairs << "\t"
        << data.mutation_count << "\n";
  }

  ofs.close();
  cout << "wrote " << bin_results.size() << " bins to " << filename << endl;
}

// Public method to orchestrate the process
void QueryBin::write_to_csv()
{
  aggregate_data();
  write_data_to_file();
}