# csegment_coverage

Compute coverage statistics for csegments (clustered segments) across multiple libraries. This command calculates the total number of sequenced base pairs (bp) that align to each csegment, where each read contributes once with its maximum intersection length across all segments in the csegment.

## Overview

The `csegment_coverage` command processes alignments from multiple libraries and computes coverage statistics for csegments. A csegment is a cluster of segments, and this command counts how many sequenced base pairs align to segments within each csegment.

The key feature is that each read contributes only once per csegment, using its maximum intersection length across all segments in that csegment. This avoids overcounting reads that map to multiple segments within the same csegment.

## Output Format

### Coverage Table (`-ofn_coverage`)

A tab-delimited table with:
- First column: `csegment` - the csegment identifier
- Subsequent columns: One column per library (using library IDs from the library table), containing total sequenced base pairs (total_bp)

The `total_bp` value represents the sum of intersection lengths for all reads that align to segments within the csegment. For each read, only its maximum intersection length (across all segments in the csegment) is included in the sum.

**Note**: The output contains raw total_bp values. These should be normalized by:
1. Mean segment length of the csegment (to get coverage: sequenced bp per bp)
2. Total reads in the library (to get CPM: coverage per million reads)

### Read Counts Table (`-ofn_read_counts`)

A tab-delimited table with:
- `lib_id`: Library identifier
- `total_reads`: Total number of reads in the alignment file for that library

## Usage

```bash
alntools csegment_coverage \
  -ifn_libraries <library_table.txt> \
  -ifn_segments <segments.txt> \
  -ifn_cseg_map <cseg_map.txt> \
  -ofn_coverage <coverage.txt> \
  -ofn_read_counts <read_counts.txt> \
  [filtering options]
```

## Parameters

### Required Input Files

- `-ifn_libraries`: Library table with `lib_id` and `aln_fn` columns specifying alignment files
- `-ifn_segments`: Segment table with `segment`, `contig`, `start`, `end` columns
- `-ifn_cseg_map`: Csegment mapping table with `segment` and `csegment` columns

### Required Output Files

- `-ofn_coverage`: Output file for coverage matrix (total_bp per csegment per library)
- `-ofn_read_counts`: Output file for read counts per library

### Filtering Parameters

- `-min_segment_length`: Minimum segment length in bp (default: 1000)
- `-clip_mode`: Alignment clipping mode (default: complete)
  - Options: `all`, `complete`, `allow_one_side_clip`, `only_one_side_clipped`, `only_two_side_clipped`, `only_clipped`, `local_align`
- `-clip_margin`: Clipping margin in bases (default: 10)
- `-min_mutations_percent`: Minimum mutations percentage for alignment filtering (default: 0.0)
- `-max_mutations_percent`: Maximum mutations percentage for alignment filtering (default: 0.1)
- `-min_alignment_length`: Minimum alignment length in read coordinates (default: 1000)
- `-max_alignment_length`: Maximum alignment length in read coordinates (default: 0, no limit)
- `-min_indel_length`: Minimum indel length to include in mutation density calculations (default: 3)

## Algorithm

For each csegment and library:

1. Identify all segments belonging to the csegment from the cseg_map
2. For each segment, find alignments that intersect the segment interval
3. Filter alignments using the specified filtering parameters
4. For each passing alignment, calculate the intersection length with the segment
5. For each read, track the maximum intersection length across all segments in the csegment
6. Sum all maximum intersection lengths to get `total_bp` for the csegment

## Normalization

The raw `total_bp` values should be normalized to get meaningful coverage metrics:

1. **Coverage** (sequenced bp per bp): `coverage = total_bp / mean_csegment_length`
   - Where `mean_csegment_length` is the average length of all segments in the csegment

2. **CPM** (coverage per million reads): `CPM = (coverage / total_reads) * 1e6`
   - This normalizes by library size to enable comparison across libraries

## Example

```bash
alntools csegment_coverage \
  -ifn_libraries libraries.txt \
  -ifn_segments segments.txt \
  -ifn_cseg_map cseg_map.txt \
  -ofn_coverage coverage.txt \
  -ofn_read_counts read_counts.txt \
  -min_segment_length 1000 \
  -clip_mode complete \
  -clip_margin 10 \
  -min_mutations_percent 0.0 \
  -max_mutations_percent 0.1 \
  -min_alignment_length 1000 \
  -min_indel_length 3
```

## Related Commands

- [cov_matrix](cov_matrix.md): Compute per-base coverage matrix for individual segments
- [coverage](coverage.md): Generate alignment coverage statistics for contigs

