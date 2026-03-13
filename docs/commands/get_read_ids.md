# get_read_ids

Assigns each read to the bin whose segment has the longest intersection with that read's alignment. Given a single ALN file and a segment table with bin identifiers, the command traverses all segments, filters alignments using the same criteria as `cov_matrix`, and outputs a table mapping every read in the ALN file to its best-matching bin (or `none` if no passing alignment intersected any segment).

## Syntax

```bash
alntools get_read_ids -ifn_aln <alignments.aln> \
                      -ifn_segments <segment_table> \
                      -ofn <output.txt> \
                      [-clip_mode <mode>] \
                      [-clip_margin <int>] \
                      [-min_mutations_percent <double>] \
                      [-max_mutations_percent <double>] \
                      [-min_alignment_length <int>] \
                      [-max_alignment_length <int>] \
                      [-min_indel_length <int>]
```

## Parameters

**Mandatory Arguments:**
* `-ifn_aln <fn>`: Input ALN file.
* `-ifn_segments <fn>`: Segment table with `contig`, `start`, `end`, `bin_id` columns.
* `-ofn <fn>`: Output tab-delimited file with `read_id` and `bin_id` columns.

**Optional Arguments - Alignment Filtering:**
* `-clip_mode <mode>`: Clipping mode for alignment filtering (default: complete). Options:
  * `all`: Allow all alignments
  * `complete`: Alignment must cover all read from start to end
  * `allow_one_side_clip`: Allow clipped on one side (start at read start or end at read end)
  * `only_one_side_clipped`: Show only alignments clipped on exactly one side
  * `only_two_side_clipped`: Show only alignments clipped on both sides
  * `only_clipped`: Show alignments clipped on one or both sides
  * `local_align`: Show only locally aligned reads (first/last alignments on same contig)
* `-clip_margin <int>`: Clipping margin in bases (default: 10). Used to determine if alignment covers read start/end.
* `-min_mutations_percent <double>`: Minimum mutations percentage (default: 0.0). Alignments with fewer mutations are excluded.
* `-max_mutations_percent <double>`: Maximum mutations percentage (default: 0.1). Alignments with more mutations are excluded.
* `-min_alignment_length <int>`: Minimum alignment length in read coordinates (default: 1000). Shorter alignments are excluded.
* `-max_alignment_length <int>`: Maximum alignment length in read coordinates (default: 0, no limit). Longer alignments are excluded when > 0.
* `-min_indel_length <int>`: Minimum indel length to include in mutation density calculations (default: 3). Smaller indels are excluded from mutation counts.

## Examples

**Basic usage with default filtering:**

```bash
alntools get_read_ids -ifn_aln output/sample.aln \
                      -ifn_segments segments.txt \
                      -ofn read_ids.txt
```

**Custom filtering for high-quality alignments:**

```bash
alntools get_read_ids -ifn_aln output/sample.aln \
                      -ifn_segments segments.txt \
                      -ofn read_ids.txt \
                      -clip_mode complete \
                      -max_mutations_percent 0.05 \
                      -min_alignment_length 2000
```

**Include all alignments without filtering:**

```bash
alntools get_read_ids -ifn_aln output/sample.aln \
                      -ifn_segments segments.txt \
                      -ofn read_ids.txt \
                      -clip_mode all \
                      -min_alignment_length 0
```

## Input File Formats

**Segment Table (`ifn_segments`):**

Tab-delimited file with header (coordinates are 1-based, inclusive):
- `contig`: Contig identifier
- `start`: Start position (1-based)
- `end`: End position (1-based, inclusive)
- `bin_id`: Bin identifier (multiple rows may share the same `bin_id`)

Example:
```
contig	start	end	bin_id
contig_1	1	5000	bin_A
contig_1	5001	10000	bin_A
contig_2	1	3500	bin_B
contig_3	2000	7000	bin_C
```

## Output File Formats

**Read ID Table (`ofn`):**

Tab-delimited file with header:
- `read_id`: Read identifier
- `bin_id`: Bin identifier of the segment with the longest intersection, or `none` if no passing alignment intersected any segment

Every read present in the ALN file appears in the output exactly once.

Example:
```
read_id	bin_id
read_001	bin_A
read_002	bin_C
read_003	none
read_004	bin_B
```

## Algorithm

1. Load the segment table
2. Load the ALN file
3. Initialize short indel counting for mutation density calculations
4. Initialize read alignment index if using `local_align` clip mode
5. For each segment, query all alignments intersecting the segment interval
6. For each passing alignment, compute the intersection length in contig coordinates
7. Update the read's assignment if this intersection is longer than the current best
8. Output all reads, writing `none` for reads with no passing intersection

When a read has multiple passing alignments to segments from different bins, the bin whose segment produces the longest single-alignment intersection wins.

## Notes

- Segment coordinates are 1-based and inclusive (standard biological convention)
- Multiple segments can share the same `bin_id`; all contribute to assignment for that bin
- Every read in the ALN file is written to the output, including reads with no intersecting alignments (`bin_id=none`)
- Alignment filtering uses the same logic as `cov_matrix`, allowing consistent filtering across commands
- Default alignment filters focus on complete, high-quality alignments (clip_mode=complete, max_mutations_percent=0.1, min_alignment_length=1000)
