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

QueryFull::QueryFull(const vector<Interval>& intervals, const AlignmentStore& store, HeightStyle height_style, int max_alignments, ClipMode clip_mode, int clip_margin, double min_mutations_percent, double max_mutations_percent, int min_alignment_length, int max_alignment_length)
    : QueryBase(intervals, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length)
    , store(store)
    , height_style(height_style)
    , max_alignments(max_alignments)
{
  // warn about performance with many intervals
  if (intervals.size() > 100) {
    cerr << "warning: QueryFull with " << intervals.size() << " intervals may be slow. Consider using fewer intervals or QueryBin for large-scale analysis." << endl;
  }
}

void QueryFull::generate_output_data()
{
  uint64_t current_alignment_index = 0;
  output_alignments.clear();
  output_mutations.clear();
  output_reads.clear();

  cout << "number of intervals: " << intervals.size() << endl;
  for (const auto& interval : intervals) {
    std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval, max_alignments);
    cout << "interval: " << interval.to_string() << endl;
    cout << "number of alignments: " << alignments.size() << endl;

    for (const auto& alignment_ref : alignments) {
      const auto& aln = alignment_ref.get();
      
      // apply alignment filtering
      if (!passes_filter(aln, store)) {
        continue;
      }
      
      string read_id = store.get_read_id(aln.read_index);
      string contig_id = store.get_contig_id(aln.contig_index);
      string cs_string = generate_cs_tag(aln, store);

      // get read length from the store
      uint32_t read_length = store.get_reads()[aln.read_index].length;

      // count mutations for this alignment
      int num_mutations = aln.get_mutation_count();

      // initialize height to 0, will be set later
      output_alignments.push_back({ current_alignment_index,
          read_id,
          static_cast<int>(read_length),
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
        // skip short indels
        if (store.is_short_indel(aln.contig_index, mutation_index)) {
          continue;
        }

        // Fetch mutation object
        const Mutation& mutation = store.get_mutation(aln.contig_index, mutation_index);

        // Position is absolute contig coordinate
        // initialize height to 0, will be set later by alignment height
        output_mutations.push_back({ current_alignment_index,
            read_id,
            contig_id,
            mutation.type,
            static_cast<int>(mutation.position + 1),
            mutation.to_string(),
            0 });
      }
      current_alignment_index++;
    }
  }

  // collect reads from alignments and calculate read heights
  collect_reads_from_alignments();
  calculate_read_heights();
  assign_alignment_heights_from_reads();
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

  // Write header with read_length and height columns
  ofs_alignments << "alignment_index\tread_id\tread_length\tcontig_id\tread_start\tread_end\tcontig_start\tcontig_end\tis_reverse\tcs_tag\tmutation_count\theight\n";

  // Write data
  for (const auto& aln_data : output_alignments) {
    ofs_alignments << aln_data.alignment_index << "\t"
                   << aln_data.read_id << "\t"
                   << aln_data.read_length << "\t"
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

  // --- Write Reads ---
  cout << "writing reads to " << ofn_prefix + "_reads.tsv" << endl;
  ofstream ofs_reads(ofn_prefix + "_reads.tsv");
  if (!ofs_reads.is_open()) {
    cerr << "error: could not open file " << ofn_prefix + "_reads.tsv" << endl;
    exit(1);
  }

  // write header
  ofs_reads << "read_id\tcontig_id\tread_length\tspan_start\tspan_end\ttotal_aligned_length\tnum_alignments\tnum_mutations\theight\tread_reversed\n";

  // write data
  for (const auto& read_data : output_reads) {
    ofs_reads << read_data.read_id << "\t"
              << read_data.contig_id << "\t"
              << read_data.read_length << "\t"
              << read_data.span_start << "\t"
              << read_data.span_end << "\t"
              << read_data.total_aligned_length << "\t"
              << read_data.num_alignments << "\t"
              << read_data.num_mutations << "\t"
              << read_data.height << "\t"
              << (read_data.read_reversed ? "true" : "false") << "\n";
  }
  ofs_reads.close();
  cout << "wrote " << output_reads.size() << " reads to " << ofn_prefix + "_reads.tsv" << endl;
}

void QueryFull::execute()
{
  // build read-to-alignments index for LOCAL_ALIGN filtering
  if (clip_mode == ClipMode::LOCAL_ALIGN) {
    // need to cast away const to call init_read_alignment_index
    const_cast<AlignmentStore&>(store).init_read_alignment_index();
  }
  
  generate_output_data();
}

void QueryFull::collect_reads_from_alignments()
{
  // group alignments by (read_id, contig_id)
  std::map<std::pair<std::string, std::string>, std::vector<const FullOutputAlignments*>> read_groups;

  for (const auto& aln : output_alignments) {
    auto key = std::make_pair(aln.read_id, aln.contig_id);
    read_groups[key].push_back(&aln);
  }

  cout << "number of unique read-contig pairs: " << read_groups.size() << endl;

  // create FullOutputReads from each group
  for (const auto& group : read_groups) {
    const std::string& read_id = group.first.first;
    const std::string& contig_id = group.first.second;
    const std::vector<const FullOutputAlignments*>& alignments = group.second;

    // calculate statistics for this read
    int span_start = alignments[0]->contig_start;
    int span_end = alignments[0]->contig_end;
    int total_aligned_length = 0;
    int num_mutations = 0;
    int read_length = alignments[0]->read_length; // should be same for all alignments

    // find the first alignment in read coordinates
    const FullOutputAlignments* first_alignment = alignments[0];
    for (const auto* aln : alignments) {
      if (aln->read_start < first_alignment->read_start) {
        first_alignment = aln;
      }
    }
    bool read_reversed = first_alignment->is_reverse;

    for (const auto* aln : alignments) {
      span_start = std::min(span_start, aln->contig_start);
      span_end = std::max(span_end, aln->contig_end);
      total_aligned_length += (aln->contig_end - aln->contig_start);
      num_mutations += aln->num_mutations;
    }

    // create read output entry with height=0 initially
    output_reads.push_back({
        read_id,
        contig_id,
        read_length,
        span_start,
        span_end,
        total_aligned_length,
        static_cast<int>(alignments.size()),
        num_mutations,
        0,
        read_reversed
    });
  }
}

void QueryFull::calculate_read_heights()
{
  if (height_style == HeightStyle::BY_COORD_LEFT) {
    calculate_heights_by_coord(true); // sort by start
  } else if (height_style == HeightStyle::BY_COORD_RIGHT) {
    calculate_heights_by_coord(false); // sort by end
  } else if (height_style == HeightStyle::BY_MUTATIONS) {
    calculate_heights_by_mutations();
  }
}

void QueryFull::assign_alignment_heights_from_reads()
{
  cout << "assigning alignment heights from reads" << endl;
  // create map from (read_id, contig_id) to read height
  std::map<std::pair<std::string, std::string>, int> read_heights;
  for (const auto& read : output_reads) {
    auto key = std::make_pair(read.read_id, read.contig_id);
    read_heights[key] = read.height;
  }

  // assign heights to alignments based on their reads
  for (auto& aln : output_alignments) {
    auto key = std::make_pair(aln.read_id, aln.contig_id);
    if (read_heights.find(key) == read_heights.end()) {
      cerr << "error: read height not found for read: " << aln.read_id << " contig: " << aln.contig_id << endl;
      exit(1);
    }
    aln.height = read_heights[key];
  }

  // assign heights to mutations based on their alignments
  std::map<uint64_t, int> alignment_heights;
  for (const auto& aln : output_alignments) {
    alignment_heights[aln.alignment_index] = aln.height;
  }


  cout << "assigning heights to mutations" << endl;
  cout << "number of mutations: " << output_mutations.size() << endl;
  for (auto& mut : output_mutations) {
    if (alignment_heights.find(mut.alignment_index) == alignment_heights.end()) {
      cerr << "error: alignment height not found for alignment: " << mut.alignment_index << endl;
      exit(1);
    }
    mut.height = alignment_heights[mut.alignment_index];
  }
  cout << "done assigning heights to " << output_mutations.size() << " mutations" << endl;
}

void QueryFull::calculate_heights_by_coord(bool sort_by_start)
{
  // group reads by contig_id
  std::map<std::string, std::vector<FullOutputReads*>> reads_by_contig;
  cout << "calculating heights by coord (" << (sort_by_start ? "left" : "right") << "), number of reads: " << output_reads.size() << endl;

  int max_coord = 0;
  for (auto& read : output_reads) {
    reads_by_contig[read.contig_id].push_back(&read);
    max_coord = std::max(max_coord, read.span_end);
  }
  cout << "Max coord: " << max_coord << endl;

  // process each contig separately
  cout << "number of contigs: " << reads_by_contig.size() << endl;
  for (auto& contig_pair : reads_by_contig) {
    auto& reads = contig_pair.second;

    // sort reads by start or end position
    if (sort_by_start) {
      std::sort(reads.begin(), reads.end(),
          [](const FullOutputReads* a, const FullOutputReads* b) {
            return a->span_start < b->span_start;
          });
    } else {
      std::sort(reads.begin(), reads.end(),
          [](const FullOutputReads* a, const FullOutputReads* b) {
            return a->span_end > b->span_end;
          });
    }

    // assign heights to avoid overlaps
    std::vector<int> height_coord; // tracks the current position at each height

    for (auto read_ptr : reads) {
      // find the lowest available height
      int height = 0;
      while (height < static_cast<int>(height_coord.size())) {
        if (sort_by_start) {
          if (read_ptr->span_start >= height_coord[height]) {
            break;
          }
         } else {
            if (read_ptr->span_end <= height_coord[height]) {
              break;
            }
        }
        height++;
      }

      // if we need a new height level
      if (height >= static_cast<int>(height_coord.size())) {
        if (sort_by_start) {
          height_coord.push_back(0);
        } else {
          height_coord.push_back(max_coord);
        }
      }

      // assign the height and update the end position
      read_ptr->height = height;
      if (sort_by_start) {
        height_coord[height] = read_ptr->span_end;
      } else {
        height_coord[height] = read_ptr->span_start;
      }
    }
  }
}

void QueryFull::calculate_heights_by_mutations()
{
  // calculate mutation density for each read
  std::vector<std::pair<int, float>> read_densities;
  cout << "calculating heights by mutations" << endl;
  cout << "number of reads: " << output_reads.size() << endl;
  for (int i = 0; i < static_cast<int>(output_reads.size()); i++) {
    const auto& read = output_reads[i];
    int aligned_length = read.total_aligned_length;
    if (aligned_length <= 0)
      aligned_length = 1; // avoid division by zero

    float density = static_cast<float>(read.num_mutations) / aligned_length;
    read_densities.push_back({ i, density });
  }

  // sort by mutation density (highest first)
  cout << "sorting read densities" << endl;
  std::sort(read_densities.begin(), read_densities.end(),
      [](const std::pair<int, float>& a, const std::pair<int, float>& b) {
        return a.second > b.second;
      });

  // group reads by contig_id for overlap prevention
  std::map<std::string, std::vector<std::vector<std::pair<int, int>>>> contig_heights;

  // assign heights in order of decreasing density while preventing overlaps
  cout << "assigning heights, number of read densities: " << read_densities.size() << endl;
  for (const auto& read_density : read_densities) {
    int read_idx = read_density.first;
    FullOutputReads& read = output_reads[read_idx];

    // get or create the heights vector for this contig
    auto& heights = contig_heights[read.contig_id];

    // find the minimum height with no overlap
    int height = 0;
    bool overlap = true;

    while (overlap) {
      // add a new height level if needed
      if (height >= static_cast<int>(heights.size())) {
        heights.push_back(std::vector<std::pair<int, int>>());
        overlap = false;
      } else {
        // check for overlaps at the current height using binary search
        const auto& intervals_at_height = heights[height];

        if (intervals_at_height.empty()) {
          // no intervals at this height yet
          overlap = false;
        } else {
          overlap = has_overlap(intervals_at_height, read.span_start, read.span_end);
        }
      }

      if (overlap) {
        height++;
      }
    }

    // assign height and add the interval to the height level
    read.height = height;

    // add the new interval and maintain sorted order
    add_sorted_interval(heights[height], read.span_start, read.span_end);
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

const std::vector<FullOutputReads>& QueryFull::get_output_reads() const
{
  return output_reads;
}

void QueryFull::set_height_style(HeightStyle style)
{
  height_style = style;
}

HeightStyle QueryFull::get_height_style() const
{
  return height_style;
}
