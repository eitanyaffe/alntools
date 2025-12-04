#include "QueryFull.h"
#include "utils.h"
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////

QueryFull::QueryFull(const vector<Interval>& intervals, const AlignmentStore& store, HeightStyle height_style, int max_alignments, ClipMode clip_mode, int clip_margin, double min_mutations_percent, double max_mutations_percent, int min_alignment_length, int max_alignment_length, int max_margin, ChunkType chunk_type)
    : QueryBase(intervals, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length)
    , store(store)
    , height_style(height_style)
    , max_alignments(max_alignments)
    , max_margin(max_margin)
    , chunk_type(chunk_type)
{
  // warn about performance with many intervals
  if (intervals.size() > 100) {
    cerr << "warning: QueryFull with " << intervals.size() << " intervals may be slow. Consider using fewer intervals or QueryBin for large-scale analysis." << endl;
  }
  
  // build map for efficient vcoord lookup
  build_contig_id_to_intervals_map();
}

////////////////////////////////////////////////////////////////////////////////
// Helper Functions
////////////////////////////////////////////////////////////////////////////////

void QueryFull::build_contig_id_to_intervals_map()
{
  contig_id_to_intervals_map.clear();
  for (const auto& interval : intervals) {
    contig_id_to_intervals_map[interval.contig].push_back(&interval);
  }
}

bool QueryFull::clip_alignment_to_interval(const Alignment& aln, const Interval& interval,
                                            int64_t& vstart, int64_t& vend,
                                            bool& clip_start, bool& clip_end) const
{
  // alignment contig coords are 0-based half-open [start, end)
  uint32_t aln_contig_start_0 = aln.contig_start;
  uint32_t aln_contig_end_0 = aln.contig_end;
  
  // check if completely outside interval (skip alignment)
  // use strict inequality to avoid warnings for boundary cases
  if (aln_contig_end_0 < interval.start || aln_contig_start_0 > interval.end) {
    return false;
  }
  
  // clip contig coordinates to interval boundaries
  uint32_t clipped_start = std::max(aln_contig_start_0, interval.start);
  uint32_t clipped_end = std::min(aln_contig_end_0, interval.end);
  
  // compute vcoords for clipped positions
  int64_t offset_start = clipped_start - interval.start;
  int64_t offset_end = clipped_end - interval.start;
  
  int64_t vstart_coord, vend_coord;
  if (interval.strand == '-') {
    // minus strand: vstart corresponds to contig end, vend to contig start
    // vcoord decreases as contig position increases
    vstart_coord = interval.vend - offset_start;
    vend_coord = interval.vend - offset_end;
    // for minus strand, vstart_coord > vend_coord, so swap
    vstart = std::min(vstart_coord, vend_coord);
    vend = std::max(vstart_coord, vend_coord);
    // clip flags are swapped for minus strand to match vcoords
    clip_start = (aln_contig_end_0 > interval.end);  // vstart maps from clipped_end
    clip_end = (aln_contig_start_0 < interval.start);  // vend maps from clipped_start
  } else {
    // plus strand: vstart corresponds to contig start
    vstart_coord = interval.vstart + offset_start;
    vend_coord = interval.vstart + offset_end;
    vstart = vstart_coord;
    vend = vend_coord;
    // clip flags match contig coordinates
    clip_start = (aln_contig_start_0 < interval.start);  // vstart maps from clipped_start
    clip_end = (aln_contig_end_0 > interval.end);  // vend maps from clipped_end
  }
  
  return true;
}

bool QueryFull::has_overlap(const std::vector<std::pair<int64_t, int64_t>>& intervals, int64_t start, int64_t end)
{
  if (intervals.empty()) {
    return false;
  }

  // Binary search to find the interval with end point just before or at the start of new interval
  auto comp = [](const std::pair<int64_t, int64_t>& interval, int64_t value) {
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

void QueryFull::add_sorted_interval(std::vector<std::pair<int64_t, int64_t>>& intervals, int64_t start, int64_t end)
{
  // Binary search to find insertion point
  auto comp = [](int64_t value, const std::pair<int64_t, int64_t>& interval) {
    return value < interval.first;
  };

  // Find the insertion point (first interval with start >= our start)
  auto it = std::upper_bound(intervals.begin(), intervals.end(), start, comp);

  // Insert the new interval
  intervals.insert(it, std::make_pair(start, end));
}

////////////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////////////

void QueryFull::write_to_csv(const std::string& odir)
{
  // --- Write Alignments ---
  cout << "writing alignments to " << odir + "/full_alignments.txt" << endl;
  ofstream ofs_alignments(odir + "/full_alignments.txt");
  if (!ofs_alignments.is_open()) {
    cerr << "error: could not open file " << odir + "/full_alignments.txt" << endl;
    exit(1);
  }

  // Write header with read_length, height, chunk_id, and strand_flipped columns
  ofs_alignments << "alignment_index\tread_id\tread_length\tcontig_id\tread_start\tread_end\tcontig_start\tcontig_end\tis_reverse\tcs_tag\tmutation_count\theight\tchunk_id\tstrand_flipped\tvstart\tvend\tclip_start\tclip_end\n";

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
                   << aln_data.height << "\t"
                   << aln_data.chunk_id << "\t"
                   << (aln_data.strand_flipped ? "true" : "false") << "\t"
                   << aln_data.vstart << "\t"
                   << aln_data.vend << "\t"
                   << (aln_data.clip_start ? "true" : "false") << "\t"
                   << (aln_data.clip_end ? "true" : "false") << "\n";
  }
  ofs_alignments.close();
  cout << "wrote " << output_alignments.size() << " alignments to " << odir + "/full_alignments.txt" << endl;

  // --- Write Mutations ---
  cout << "writing mutations to " << odir + "/full_mutations.txt" << endl;
  ofstream ofs_mutations(odir + "/full_mutations.txt");
  if (!ofs_mutations.is_open()) {
    cerr << "error: could not open file " << odir + "/full_mutations.txt" << endl;
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
  cout << "wrote " << output_mutations.size() << " mutations to " << odir + "/full_mutations.txt" << endl;

  // --- Write Chunks ---
  cout << "writing chunks to " << odir + "/full_chunks.txt" << endl;
  ofstream ofs_chunks(odir + "/full_chunks.txt");
  if (!ofs_chunks.is_open()) {
    cerr << "error: could not open file " << odir + "/full_chunks.txt" << endl;
    exit(1);
  }

  // write header
  ofs_chunks << "chunk_id\tread_id\tread_length\ttotal_aligned_length\tnum_alignments\tnum_mutations\theight\tread_reversed\tvstart\tvend\n";

  // write data
  for (const auto& chunk_data : output_chunks) {
    ofs_chunks << chunk_data.chunk_id << "\t"
               << chunk_data.read_id << "\t"
               << chunk_data.read_length << "\t"
               << chunk_data.total_aligned_length << "\t"
               << chunk_data.num_alignments << "\t"
               << chunk_data.num_mutations << "\t"
               << chunk_data.height << "\t"
               << (chunk_data.read_reversed ? "true" : "false") << "\t"
               << chunk_data.vstart << "\t"
               << chunk_data.vend << "\n";
  }
  ofs_chunks.close();
  cout << "wrote " << output_chunks.size() << " chunks to " << odir + "/full_chunks.txt" << endl;
}

const std::vector<FullOutputAlignments>& QueryFull::get_output_alignments() const
{
  return output_alignments;
}

const std::vector<FullOutputMutations>& QueryFull::get_output_mutations() const
{
  return output_mutations;
}

const std::vector<FullOutputChunk>& QueryFull::get_output_chunks() const
{
  return output_chunks;
}

void QueryFull::set_height_style(HeightStyle style)
{
  height_style = style;
}

HeightStyle QueryFull::get_height_style() const
{
  return height_style;
  }

////////////////////////////////////////////////////////////////////////////////
// Main Pipeline Functions
////////////////////////////////////////////////////////////////////////////////

void QueryFull::collect_alignments_from_intervals()
{
  uint64_t current_alignment_index = 0;
  output_alignments.clear();
  output_mutations.clear();

  cout << "number of intervals: " << intervals.size() << endl;
  
  // calculate total interval length for proportional max_alignments distribution
  uint64_t total_length = 0;
  for (const auto& interval : intervals) {
    total_length += (interval.end - interval.start);
  }
  
  for (const auto& interval : intervals) {
    // calculate interval-specific max_alignments based on relative length
    int interval_max_alignments = max_alignments;
    if (max_alignments > 0 && total_length > 0) {
      uint64_t interval_length = interval.end - interval.start;
      double weight = static_cast<double>(interval_length) / static_cast<double>(total_length);
      interval_max_alignments = static_cast<int>(max_alignments * weight);
      // ensure at least 1 if max_alignments > 0 and interval has length
      if (interval_max_alignments == 0 && interval_length > 0) {
        interval_max_alignments = 1;
      }
    }
    
    std::vector<std::reference_wrapper<const Alignment>> alignments = 
      store.get_alignments_intersecting_interval(interval, interval_max_alignments);
    cout << "interval: " << interval.to_string() << " (max_alignments: " << interval_max_alignments << ")" << endl;
    cout << "number of alignments: " << alignments.size() << endl;

    for (const auto& alignment_ref : alignments) {
      const auto& aln = alignment_ref.get();
      
      string read_id = store.get_read_id(aln.read_index);
      
      // apply alignment filtering
      if (!passes_filter(aln, store)) {
        continue;
      }
      
      string contig_id = store.get_contig_id(aln.contig_index);
      string cs_string = generate_cs_tag(aln, store);

      // get read length from the store
      uint32_t read_length = store.get_reads()[aln.read_index].length;

      // count mutations for this alignment
      int num_mutations = aln.get_mutation_count();

      // clip alignment to interval and compute vcoords
      int64_t vstart, vend;
      bool clip_start, clip_end;
      
      if (!clip_alignment_to_interval(aln, interval, vstart, vend, clip_start, clip_end)) {
        // alignment completely outside interval, skip it
        continue;
      }
      
      // filter based on clipped alignment length
      uint32_t clipped_start = std::max(aln.contig_start, interval.start);
      uint32_t clipped_end = std::min(aln.contig_end, interval.end);
      uint32_t clipped_length = clipped_end - clipped_start;
      
      if (min_alignment_length > 0 && clipped_length < static_cast<uint32_t>(min_alignment_length)) {
        continue;
      }

      // initialize height to 0, will be set later
      // strand_flipped will be set later when building chunks
      // convert to 1-based coordinates for output
      output_alignments.push_back({ current_alignment_index,
          read_id,
          static_cast<int>(read_length),
          contig_id,
          static_cast<int>(aln.read_start + 1),
          static_cast<int>(aln.read_end),
          static_cast<int>(aln.contig_start + 1),
          static_cast<int>(aln.contig_end),
          aln.is_reverse,
          cs_string,
          num_mutations,
          0,
          "",
          false,  // strand_flipped will be set in build_chunk_objects
          vstart,
          vend,
          clip_start,
          clip_end });

      for (uint32_t mutation_index : aln.mutations) { // Iterate indices
        // skip short indels
        if (should_skip_short_indel(store, aln.contig_index, mutation_index)) {
          continue;
}

        // Fetch mutation object
        const Mutation& mutation = store.get_mutation(aln.contig_index, mutation_index);
    
        // Position is absolute contig coordinate (1-based for output)
        uint32_t position_1based = mutation.position + 1;
        

        // initialize height to 0, will be set later by alignment height
        output_mutations.push_back({ current_alignment_index,
            read_id,
            contig_id,
            mutation.type,
            static_cast<int>(position_1based),
            mutation.to_string(),
            0 });
            
      }
      current_alignment_index++;
    }
  }
}

std::vector<std::vector<FullOutputAlignments*>> QueryFull::group_alignments_into_chunks(
  const std::vector<FullOutputAlignments*>& alignments)
{
    std::vector<std::vector<FullOutputAlignments*>> chunks;
    
    if (chunk_type == ChunkType::READ) {
      // entire read is one chunk
      if (!alignments.empty()) {
        chunks.push_back(alignments);
      }
    return chunks;
  }
  
  if (chunk_type == ChunkType::ALIGNMENT) {
      // each alignment is its own chunk
      for (auto* aln : alignments) {
        chunks.push_back({aln});
      }
    return chunks;
  }
  
      // sort alignments by read_start for overlap and gap-based chunking
  std::vector<FullOutputAlignments*> sorted_alignments = alignments;
  std::sort(sorted_alignments.begin(), sorted_alignments.end(),
          [](const FullOutputAlignments* a, const FullOutputAlignments* b) {
            return a->read_start < b->read_start;
          });

  if (sorted_alignments.empty()) {
    return chunks;
  }

      if (chunk_type == ChunkType::BREAK_ON_OVERLAP) {
        // start new chunk if next alignment overlaps: next.start < (current_chunk_end - max_margin)
    chunks.push_back({sorted_alignments[0]});
    int current_chunk_end = sorted_alignments[0]->read_end;
          
    for (size_t i = 1; i < sorted_alignments.size(); i++) {
            // check if next alignment overlaps
      if (sorted_alignments[i]->read_start < (current_chunk_end - max_margin)) {
              // start new chunk due to overlap
        chunks.push_back({sorted_alignments[i]});
        current_chunk_end = sorted_alignments[i]->read_end;
            } else {
              // add to current chunk and update end coordinate
        chunks.back().push_back(sorted_alignments[i]);
        current_chunk_end = std::max(current_chunk_end, sorted_alignments[i]->read_end);
            }
          }
  } else { // ChunkType::BREAK_ON_GAP
        // break if next alignment has a large gap
    chunks.push_back({sorted_alignments[0]});
          
    for (size_t i = 1; i < sorted_alignments.size(); i++) {
      int gap = sorted_alignments[i]->read_start - sorted_alignments[i-1]->read_end;
            if (gap <= max_margin) {
              // add to current chunk
        chunks.back().push_back(sorted_alignments[i]);
            } else {
              // start new chunk
        chunks.push_back({sorted_alignments[i]});
            }
          }
        }
  
  return chunks;
    }

void QueryFull::build_chunk_objects(const std::string& read_id, 
  const std::vector<std::vector<FullOutputAlignments*>>& chunks)
{
    for (size_t chunk_idx = 0; chunk_idx < chunks.size(); chunk_idx++) {
      const auto& chunk_alignments = chunks[chunk_idx];
      
      // generate chunk_id
      std::string chunk_id = read_id + "_chunk_" + std::to_string(chunk_idx);
      
      // calculate chunk statistics using vcoords
      int64_t chunk_vstart = INT64_MAX;
      int64_t chunk_vend = INT64_MIN;
      int total_aligned_length = 0;
      int num_mutations = 0;
      int read_length = chunk_alignments[0]->read_length;
      
      // find the first alignment in read coordinates to determine read direction
      const FullOutputAlignments* first_alignment = chunk_alignments[0];
      for (const auto* aln : chunk_alignments) {
        if (aln->read_start < first_alignment->read_start) {
          first_alignment = aln;
        }
      }
      bool read_reversed = first_alignment->is_reverse;
      std::string first_alignment_contig = first_alignment->contig_id;
      bool first_alignment_strand = first_alignment->is_reverse;

      for (auto* aln : chunk_alignments) {
        // set strand_flipped: true if same contig as first AND different strand
        aln->strand_flipped = (aln->contig_id == first_alignment_contig) && 
                              (aln->is_reverse != first_alignment_strand);
        // use vcoords already computed during clip_alignment_to_interval
        // these are already clipped to interval boundaries and correctly handle strand
        int64_t aln_vstart = aln->vstart;
        int64_t aln_vend = aln->vend;
        
        if (aln_vstart >= 0 && aln_vend >= 0) {
          int64_t aln_min = std::min(aln_vstart, aln_vend);
          int64_t aln_max = std::max(aln_vstart, aln_vend);
          // take min/max of both vcoords (handles minus strand where vstart > vend)
          chunk_vstart = std::min(chunk_vstart, aln_min);
          chunk_vend = std::max(chunk_vend, aln_max);
        }
        
        total_aligned_length += (aln->contig_end - aln->contig_start);
        num_mutations += aln->num_mutations;
      }

      // handle case where no valid vcoords found
      if (chunk_vstart == INT64_MAX) chunk_vstart = 0;
      if (chunk_vend == INT64_MIN) chunk_vend = 0;

      // create chunk output entry with vcoords
      output_chunks.push_back({
          chunk_id,
          read_id,
          read_length,
          total_aligned_length,
          static_cast<int>(chunk_alignments.size()),
          num_mutations,
          0, // height will be set later
          read_reversed,
          chunk_vstart,
          chunk_vend
      });

      // assign chunk_id to alignments and update mapping
      int chunk_index = static_cast<int>(output_chunks.size()) - 1;
      for (auto* aln : chunk_alignments) {
        aln->chunk_id = chunk_id;
        alignment_to_chunk_index[aln->alignment_index] = chunk_index;
      }
    }
}

void QueryFull::collect_chunks_from_alignments()
{
  // group alignments by read_id only (not read_id + contig_id)
  std::map<std::string, std::vector<FullOutputAlignments*>> read_groups;

  for (auto& aln : output_alignments) {
    read_groups[aln.read_id].push_back(&aln);
  }

  cout << "creating chunks from " << read_groups.size() << " unique reads using chunk_type ";
  switch (chunk_type) {
    case ChunkType::READ: cout << "READ"; break;
    case ChunkType::ALIGNMENT: cout << "ALIGNMENT"; break;
    case ChunkType::BREAK_ON_OVERLAP: cout << "BREAK_ON_OVERLAP"; break;
    case ChunkType::BREAK_ON_GAP: cout << "BREAK_ON_GAP"; break;
  }
  cout << endl;

  // create chunks from each read group
  alignment_to_chunk_index.resize(output_alignments.size(), -1);

  for (auto& group : read_groups) {
    const std::string& read_id = group.first;
    auto& alignments = group.second;

    std::vector<std::vector<FullOutputAlignments*>> chunks = group_alignments_into_chunks(alignments);
    build_chunk_objects(read_id, chunks);
  }

  cout << "created " << output_chunks.size() << " chunks" << endl;
}

void QueryFull::calculate_chunk_heights()
{
  if (height_style == HeightStyle::BY_COORD_LEFT) {
    calculate_chunk_heights_by_coord(true); // sort by start
  } else if (height_style == HeightStyle::BY_COORD_RIGHT) {
    calculate_chunk_heights_by_coord(false); // sort by end
  } else if (height_style == HeightStyle::BY_MUTATIONS) {
    calculate_chunk_heights_by_mutations();
  }
}

void QueryFull::calculate_chunk_heights_by_coord(bool sort_by_start)
{
  // process all chunks globally using vcoords (which span contigs)
  cout << "calculating chunk heights by coord (" << (sort_by_start ? "left" : "right") << "), number of chunks: " << output_chunks.size() << endl;

  // collect all chunk pointers
  std::vector<FullOutputChunk*> chunks;
  int64_t max_coord = 0;
  for (auto& chunk : output_chunks) {
    chunks.push_back(&chunk);
    max_coord = std::max(max_coord, chunk.vend);
  }

  // sort chunks by vstart or vend position globally
  if (sort_by_start) {
    std::sort(chunks.begin(), chunks.end(),
        [](const FullOutputChunk* a, const FullOutputChunk* b) {
          return a->vstart < b->vstart;
        });
  } else {
    std::sort(chunks.begin(), chunks.end(),
        [](const FullOutputChunk* a, const FullOutputChunk* b) {
          return a->vend > b->vend;
        });
  }

  // assign heights to avoid overlaps globally using vcoords
  std::vector<int64_t> height_coord; // tracks the current position at each height

  for (auto chunk_ptr : chunks) {
    // find the lowest available height
    int height = 0;
    while (height < static_cast<int>(height_coord.size())) {
      if (sort_by_start) {
        if (chunk_ptr->vstart >= height_coord[height]) {
          break;
        }
       } else {
          if (chunk_ptr->vend <= height_coord[height]) {
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
    chunk_ptr->height = height;
    if (sort_by_start) {
      height_coord[height] = chunk_ptr->vend;
    } else {
      height_coord[height] = chunk_ptr->vstart;
    }
  }
}

void QueryFull::calculate_chunk_heights_by_mutations()
{
  // calculate mutation density and multiple alignment status for each chunk
  std::vector<std::tuple<int, bool, float>> chunk_data;
  cout << "calculating chunk heights by mutations" << endl;
  cout << "number of chunks: " << output_chunks.size() << endl;
  for (int i = 0; i < static_cast<int>(output_chunks.size()); i++) {
    const auto& chunk = output_chunks[i];
    int aligned_length = chunk.total_aligned_length;
    if (aligned_length <= 0)
      aligned_length = 1; // avoid division by zero

    float density = static_cast<float>(chunk.num_mutations) / aligned_length;
    bool has_multiple_alignments = chunk.num_alignments > 1;
    chunk_data.push_back(std::make_tuple(i, has_multiple_alignments, density));
  }

  // sort by has_multiple_alignments (single alignments first), then by mutation density (highest first)
  cout << "sorting chunks by multiple alignments and mutation density" << endl;
  std::sort(chunk_data.begin(), chunk_data.end(),
      [](const std::tuple<int, bool, float>& a, const std::tuple<int, bool, float>& b) {
        // primary sort: has_multiple_alignments (false first, true last)
        if (std::get<1>(a) != std::get<1>(b)) {
          return std::get<1>(a) < std::get<1>(b);
        }
        // secondary sort: mutation density (highest first)
        return std::get<2>(a) > std::get<2>(b);
      });

  // use global heights for overlap prevention (vcoords span contigs)
  std::vector<std::vector<std::pair<int64_t, int64_t>>> global_heights;

  // assign heights in sorted order while preventing overlaps globally
  cout << "assigning heights, number of chunks: " << chunk_data.size() << endl;
  for (const auto& chunk_entry : chunk_data) {
    int chunk_idx = std::get<0>(chunk_entry);
    FullOutputChunk& chunk = output_chunks[chunk_idx];

    // find the minimum height with no overlap globally
    int height = 0;
    bool overlap = true;

    while (overlap) {
      // add a new height level if needed
      if (height >= static_cast<int>(global_heights.size())) {
        global_heights.push_back(std::vector<std::pair<int64_t, int64_t>>());
        overlap = false;
      } else {
        // check for overlaps at the current height using binary search
        const auto& intervals_at_height = global_heights[height];

        if (intervals_at_height.empty()) {
          // no intervals at this height yet
          overlap = false;
        } else {
          overlap = has_overlap(intervals_at_height, chunk.vstart, chunk.vend);
        }
      }

      if (overlap) {
        height++;
      }
    }

    // assign height and add the interval to the height level
    chunk.height = height;

    // add the new interval and maintain sorted order
    add_sorted_interval(global_heights[height], chunk.vstart, chunk.vend);
  }
}

void QueryFull::assign_alignment_heights_from_chunks()
{
  cout << "assigning alignment heights from chunks" << endl;
  
  // assign heights to alignments based on their chunks
  for (auto& aln : output_alignments) {
    if (alignment_to_chunk_index[aln.alignment_index] == -1) {
      cerr << "error: chunk index not found for alignment: " << aln.alignment_index << endl;
      exit(1);
    }
    int chunk_index = alignment_to_chunk_index[aln.alignment_index];
    aln.height = output_chunks[chunk_index].height;
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

////////////////////////////////////////////////////////////////////////////////
// Verification
////////////////////////////////////////////////////////////////////////////////

void QueryFull::verify_chunks_and_alignments() const
{
  cout << "verifying chunks and alignments consistency..." << endl;
  
  // build map from chunk_id to chunk
  std::map<std::string, const FullOutputChunk*> chunk_map;
  for (const auto& chunk : output_chunks) {
    if (chunk_map.find(chunk.chunk_id) != chunk_map.end()) {
      cerr << "WARNING: duplicate chunk_id found: " << chunk.chunk_id << endl;
    }
    chunk_map[chunk.chunk_id] = &chunk;
  }
  
  // build map from chunk_id to list of alignments
  std::map<std::string, std::vector<const FullOutputAlignments*>> chunk_to_alignments;
  for (const auto& aln : output_alignments) {
    if (aln.chunk_id.empty()) {
      cerr << "WARNING: alignment " << aln.alignment_index 
           << " (read_id: " << aln.read_id 
           << ") has empty chunk_id" << endl;
    } else {
      chunk_to_alignments[aln.chunk_id].push_back(&aln);
    }
  }
  
  // verify every chunk has at least one alignment
  bool has_warnings = false;
  for (const auto& chunk : output_chunks) {
    auto it = chunk_to_alignments.find(chunk.chunk_id);
    if (it == chunk_to_alignments.end() || it->second.empty()) {
      cerr << "WARNING: chunk " << chunk.chunk_id 
           << " (read_id: " << chunk.read_id 
           << ", vstart: " << chunk.vstart 
           << ", vend: " << chunk.vend 
           << ", num_alignments: " << chunk.num_alignments 
           << ") has no associated alignments" << endl;
      has_warnings = true;
    } else {
      // verify chunk coordinates span the alignments
      int64_t min_aln_vstart = INT64_MAX;
      int64_t max_aln_vend = INT64_MIN;
      bool found_valid_vcoords = false;
      
      for (const auto* aln : it->second) {
        // handle both vstart and vend (for minus strand, vstart > vend)
        int64_t aln_min = std::min(aln->vstart, aln->vend);
        int64_t aln_max = std::max(aln->vstart, aln->vend);
        
        if (aln->vstart >= 0 && aln->vend >= 0) {
          found_valid_vcoords = true;
          min_aln_vstart = std::min(min_aln_vstart, aln_min);
          max_aln_vend = std::max(max_aln_vend, aln_max);
        }
      }
      
      if (found_valid_vcoords) {
        // chunk vcoords should span alignment vcoords (with some tolerance for rounding)
        int64_t chunk_min = std::min(chunk.vstart, chunk.vend);
        int64_t chunk_max = std::max(chunk.vstart, chunk.vend);
        
        if (chunk_min > min_aln_vstart || chunk_max < max_aln_vend) {
          cerr << "WARNING: chunk " << chunk.chunk_id 
               << " (read_id: " << chunk.read_id 
               << ") coordinates do not span its alignments:" << endl;
          cerr << "  chunk vcoords: [" << chunk.vstart << ", " << chunk.vend << "]" << endl;
          cerr << "  chunk min/max: [" << chunk_min << ", " << chunk_max << "]" << endl;
          cerr << "  alignment vcoords range: [" << min_aln_vstart << ", " << max_aln_vend << "]" << endl;
          cerr << "  number of alignments: " << it->second.size() << endl;
          cerr << "  alignment details:" << endl;
          for (const auto* aln : it->second) {
            cerr << "    alignment_index: " << aln->alignment_index 
                 << ", vstart: " << aln->vstart 
                 << ", vend: " << aln->vend 
                 << ", contig: " << aln->contig_id << endl;
          }
          has_warnings = true;
        }
      }
    }
  }
  
  // verify every alignment has a valid chunk
  for (const auto& aln : output_alignments) {
    if (!aln.chunk_id.empty() && chunk_map.find(aln.chunk_id) == chunk_map.end()) {
      cerr << "WARNING: alignment " << aln.alignment_index 
           << " (read_id: " << aln.read_id 
           << ", chunk_id: " << aln.chunk_id 
           << ") references non-existent chunk" << endl;
      has_warnings = true;
    }
  }
  
  if (has_warnings) {
    cerr << "WARNING: chunks and alignments have inconsistencies (continuing anyway)" << endl;
  } else {
    cout << "verification passed: " << output_chunks.size() << " chunks and " 
         << output_alignments.size() << " alignments are consistent" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Execute
////////////////////////////////////////////////////////////////////////////////

void QueryFull::execute()
{
  // build read-to-alignments index for LOCAL_ALIGN filtering
  init_local_align_if_needed(const_cast<AlignmentStore&>(store), clip_mode);
  
  // collect alignments and mutations from intervals
  collect_alignments_from_intervals();
  
  // group alignments into chunks based on chunk_type (READ, ALIGNMENT, BREAK_ON_OVERLAP, BREAK_ON_GAP)
  collect_chunks_from_alignments();
  
  // assign heights to chunks based on height_style (BY_COORD_LEFT, BY_COORD_RIGHT, BY_MUTATIONS)
  calculate_chunk_heights();
  
  // propagate chunk heights to their constituent alignments
  assign_alignment_heights_from_chunks();
  
  // verify chunks and alignments are properly matched
  verify_chunks_and_alignments();
}

