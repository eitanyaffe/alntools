#include "QueryFull.h"
#include "utils.h"
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

QueryFull::QueryFull(const vector<Interval>& intervals, const AlignmentStore& store, HeightStyle height_style)
    : intervals(intervals)
    , store(store)
    , height_style(height_style)
{
}

void QueryFull::generate_output_data()
{
  uint64_t current_alignment_index = 0;
  output_alignments.clear();
  output_mutations.clear();

  cout << "number of intervals: " << intervals.size() << endl;
  for (const auto& interval : intervals) {
    std::vector<std::reference_wrapper<const Alignment>> alignments =
      store.get_alignments_in_interval(interval);
    cout << "interval: " << interval.to_string() << endl;
    cout << "number of alignments: " << alignments.size() << endl;
    for (const auto& alignment_ref : alignments) {
      const auto& aln = alignment_ref.get();
      string read_id = store.get_read_id(aln.read_index);
      string contig_id = store.get_contig_id(aln.contig_index);
      string cs_string = generate_cs_tag(aln, store);

      // Count mutations for this alignment
      int num_mutations = aln.mutations.size();

      // initialize height to 0, will be set later
      output_alignments.push_back({ current_alignment_index,
          read_id,
          contig_id,
          static_cast<int>(aln.read_start),
          static_cast<int>(aln.read_end),
          static_cast<int>(aln.contig_start),
          static_cast<int>(aln.contig_end),
          aln.is_reverse,
          cs_string,
          num_mutations,
          0 });

      for (uint32_t mutation_index : aln.mutations) { // Iterate indices
        // Fetch mutation object
        const Mutation& mutation = store.get_mutation(aln.contig_index, mutation_index);

        // Position is absolute contig coordinate
        // initialize height to 0, will be set later by alignment height
        output_mutations.push_back({ current_alignment_index,
            read_id,
            contig_id,
            mutation.type,
            static_cast<int>(mutation.position),
            mutation.to_string(),
            0 });
      }
      current_alignment_index++;
    }
  }

  // calculate heights after all alignments are collected
  calculate_heights();
}

void QueryFull::write_to_csv(const std::string& ofn_prefix)
{
  // --- Write Alignments ---
  cout << "writing alignments to " << ofn_prefix + "_alignments.tsv" << endl;
  ofstream ofs_alignments(ofn_prefix + "_alignments.tsv");
  if (!ofs_alignments.is_open()) {
    cerr << "error: could not open file " << ofn_prefix + "_alignments.tsv" << endl;
    exit(1);
  }

  // Write header with added height column
  ofs_alignments << "alignment_index\tread_id\tcontig_id\tread_start\tread_end\tcontig_start\tcontig_end\tis_reverse\tcs_tag\tmutation_count\theight\n";

  // Write data
  for (const auto& aln_data : output_alignments) {
    ofs_alignments << aln_data.alignment_index << "\t"
                   << aln_data.read_id << "\t"
                   << aln_data.contig_id << "\t"
                   << aln_data.read_start << "\t"
                   << aln_data.read_end << "\t"
                   << aln_data.contig_start << "\t"
                   << aln_data.contig_end << "\t"
                   << (aln_data.is_reverse ? "true" : "false") << "\t"
                   << aln_data.cs_tag << "\t"
                   << aln_data.num_mutations << "\t"
                   << aln_data.height << "\n";
  }
  ofs_alignments.close();
  cout << "wrote " << output_alignments.size() << " alignments to " << ofn_prefix + "_alignments.tsv" << endl;

  // --- Write Mutations ---
  cout << "writing mutations to " << ofn_prefix + "_mutations.tsv" << endl;
  ofstream ofs_mutations(ofn_prefix + "_mutations.tsv");
  if (!ofs_mutations.is_open()) {
    cerr << "error: could not open file " << ofn_prefix + "_mutations.tsv" << endl;
    exit(1);
  }

  // Write header with added height column
  ofs_mutations << "alignment_index\tread_id\tcontig_id\tmutation_type\tmutation_position\tmutation_desc\theight\n";

  // Write data
  for (const auto& mut_data : output_mutations) {
    ofs_mutations << mut_data.alignment_index << "\t"
                  << mut_data.read_id << "\t"
                  << mut_data.contig_id << "\t"
                  << mut_data.type << "\t"
                  << mut_data.position << "\t"
                  << mut_data.desc << "\t"
                  << mut_data.height << "\n";
  }
  ofs_mutations.close();
  cout << "wrote " << output_mutations.size() << " mutations to " << ofn_prefix + "_mutations.tsv" << endl;
}

void QueryFull::execute()
{
  generate_output_data();
}

void QueryFull::calculate_heights()
{
  if (height_style == HeightStyle::BY_COORD) {
    calculate_heights_by_coord();
  } else if (height_style == HeightStyle::BY_MUTATIONS) {
    calculate_heights_by_mutations();
  }

  // update heights for mutations based on their alignment heights
  std::map<uint64_t, int> alignment_heights;
  for (const auto& aln : output_alignments) {
    alignment_heights[aln.alignment_index] = aln.height;
  }

  for (auto& mut : output_mutations) {
    if (alignment_heights.find(mut.alignment_index) != alignment_heights.end()) {
      mut.height = alignment_heights[mut.alignment_index];
    }
  }
}

void QueryFull::calculate_heights_by_coord()
{
  // Group alignments by contig_id
  std::map<std::string, std::vector<FullOutputAlignments*>> alignments_by_contig;
  cout << "calculating heights by coord, number of alignments: " << output_alignments.size() << endl;

  for (auto& aln : output_alignments) {
    alignments_by_contig[aln.contig_id].push_back(&aln);
  }

  // Process each contig separately
  cout << "number of contigs: " << alignments_by_contig.size() << endl;
  for (auto& contig_pair : alignments_by_contig) {
    auto& alignments = contig_pair.second;

    // Sort alignments by start position
    std::sort(alignments.begin(), alignments.end(),
        [](const FullOutputAlignments* a, const FullOutputAlignments* b) {
          return a->contig_start < b->contig_start;
        });

    // Assign heights to avoid overlaps
    std::vector<int> height_ends; // tracks the end position at each height

    for (auto aln_ptr : alignments) {
      // Find the lowest available height
      int height = 0;
      while (height < static_cast<int>(height_ends.size())) {
        if (aln_ptr->contig_start >= height_ends[height]) {
          break;
        }
        height++;
      }

      // If we need a new height level
      if (height >= static_cast<int>(height_ends.size())) {
        height_ends.push_back(0);
      }

      // Assign the height and update the end position
      aln_ptr->height = height;
      height_ends[height] = aln_ptr->contig_end;
    }
  }
}

void QueryFull::calculate_heights_by_mutations()
{
  // Calculate mutation density for each alignment
  std::vector<std::pair<int, float>> alignment_densities;
  cout << "calculating heights by mutations" << endl;
  // Calculate densities - we can now use the num_mutations field directly
  cout << "number of alignments: " << output_alignments.size() << endl;
  for (int i = 0; i < static_cast<int>(output_alignments.size()); i++) {
    const auto& aln = output_alignments[i];
    int aln_length = aln.contig_end - aln.contig_start;
    if (aln_length <= 0)
      aln_length = 1; // avoid division by zero

    float density = static_cast<float>(aln.num_mutations) / aln_length;
    alignment_densities.push_back({ i, density });
  }

  // Sort by mutation density (highest first)
  cout << "sorting alignment densities" << endl;
  std::sort(alignment_densities.begin(), alignment_densities.end(),
      [](const std::pair<int, float>& a, const std::pair<int, float>& b) {
        return a.second > b.second;
      });

  // Group alignments by contig_id for overlap prevention
  std::map<std::string, std::vector<std::vector<std::pair<int, int>>>> contig_heights;

  // Assign heights in order of decreasing density while preventing overlaps
  cout << "assigning heights, number of mutation densities: " << alignment_densities.size() << endl;
  for (const auto& aln_density : alignment_densities) {
    int aln_idx = aln_density.first;
    FullOutputAlignments& aln = output_alignments[aln_idx];

    // Get or create the heights vector for this contig
    auto& heights = contig_heights[aln.contig_id];

    // Find the minimum height with no overlap
    int height = 0;
    bool overlap = true;

    while (overlap) {
      // Add a new height level if needed
      if (height >= static_cast<int>(heights.size())) {
        heights.push_back(std::vector<std::pair<int, int>>());
        overlap = false;
      } else {
        // Check for overlaps at the current height using binary search
        const auto& intervals_at_height = heights[height];

        if (intervals_at_height.empty()) {
          // No intervals at this height yet
          overlap = false;
        } else {
          overlap = has_overlap(intervals_at_height, aln.contig_start, aln.contig_end);
        }
      }

      if (overlap) {
        height++;
      }
    }

    // Assign height and add the interval to the height level
    aln.height = height;

    // Add the new interval and maintain sorted order
    add_sorted_interval(heights[height], aln.contig_start, aln.contig_end);
  }
}

// helper to check if a new interval overlaps with any existing intervals using binary search
bool QueryFull::has_overlap(const std::vector<std::pair<int, int>>& intervals, int start, int end)
{
  if (intervals.empty()) {
    return false;
  }

  // Binary search to find the interval with end point just before or at the start of new interval
  auto comp = [](const std::pair<int, int>& interval, int value) {
    return interval.second < value;
  };

  // Find the first interval whose end >= start (might overlap)
  auto it = std::lower_bound(intervals.begin(), intervals.end(), start, comp);

  // If we found an interval and it overlaps
  if (it != intervals.end()) {
    // Check if the found interval overlaps
    if (it->first < end) {
      return true;
    }
  }

  // Check the previous interval too (if exists) since it might extend into our new interval
  if (it != intervals.begin()) {
    auto prev = it - 1;
    if (prev->second > start) {
      return true;
    }
  }

  return false;
}

// helper to add an interval to a sorted vector of intervals
void QueryFull::add_sorted_interval(std::vector<std::pair<int, int>>& intervals, int start, int end)
{
  // Binary search to find insertion point
  auto comp = [](int value, const std::pair<int, int>& interval) {
    return value < interval.first;
  };

  // Find the insertion point (first interval with start >= our start)
  auto it = std::upper_bound(intervals.begin(), intervals.end(), start, comp);

  // Insert the new interval
  intervals.insert(it, std::make_pair(start, end));
}

const std::vector<FullOutputAlignments>& QueryFull::get_output_alignments() const
{
  return output_alignments;
}

const std::vector<FullOutputMutations>& QueryFull::get_output_mutations() const
{
  return output_mutations;
}

void QueryFull::set_height_style(HeightStyle style)
{
  height_style = style;
}

HeightStyle QueryFull::get_height_style() const
{
  return height_style;
}
