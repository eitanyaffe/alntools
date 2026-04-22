# seg_matrix

Generates segment adjacency and reach matrices by analyzing read alignments. The adjacency matrix counts reads supporting direct adjacency between oriented segment sides when their estimated read-space gap is within the distance threshold, while the reach matrix counts reads supporting connections between segment sides without the adjacency constraint.

## Syntax

```bash
alntools seg_matrix -ifn_aln <input.aln> \
                    -ifn_segments <segment_table> \
                    -odir <output_dir> \
                    [options...]
```

## Parameters

**Mandatory Arguments:**
* `-ifn_aln <fn>`: Input ALN file
* `-ifn_segments <fn>`: Segment table with `segment`, `contig`, `start`, `end` columns (tab-delimited, 1-based coordinates).
* `-odir <dir>`: Output directory for result files.

**Optional Arguments:**
* `-max_mutation_percent <double>`: Maximum mutation percentage (default: 1.0). Used for (1) discarding whole alignments whose span exceeds this value, and (2) requiring covered sides and association endpoints to be at or below this value. **Definition:** for a contig interval `[start, end)`, the value is `(number of mutation events with position in that interval) / (end - start) * 100`. Each indel counts as one event regardless of its length (so one large insertion is one event). Mutations with indel length below `min_indel_length` are excluded from this count.
* `-max_adjacency_distance <int>`: Maximum distance in bp for adjacency (default: 200). Distance is measured in estimated read space.
* `-max_margin <int>`: Maximum overlap between consecutive alignments in bp (default: 10). Alignments overlapping more than this threshold are skipped.
* `-min_indel_length <int>`: Minimum indel length to include in mutation density calculations (default: 3). Shorter indels are excluded.
* `-max_side_indel_bp <int>`: Maximum total indel bp within a side region to accept coverage (default: 100). If the combined length of all deletions overlapping the side region plus insertions anchored within it (each with length ≥ `min_indel_length`) exceeds this threshold, the side is not marked as covered.
* `-side_length <int>`: Length in bp of each segment side used for coverage analysis (default: 1000).
* `-side_margin <int>`: Margin in bp between the segment edge and the start of each side region (default: 200).
* `-output_read_details <bool>`: If true, write `read_associations_adjacency.txt` and `read_associations_reach.txt` under `odir` (default: false).
* `-ifn_segment_clusters <fn>`: Optional segment-cluster table with `segment`, `csegment`, `strand` columns. When set, every segment in `ifn_segments` must appear in this file; the command also writes `cluster_adjacency_matrix.txt`, `cluster_reach_matrix.txt`, and `cluster_summary.txt`.

## Input File Formats

**Segment Table (`ifn_segments`):**

Tab-delimited file with header (coordinates are 1-based, inclusive):
* `segment`: Unique segment identifier
* `contig`: Contig identifier
* `start`: Start position (1-based)
* `end`: End position (1-based, inclusive)

Example:
```
segment	contig	start	end
seg_001	contig_1	1	5000
seg_002	contig_1	5001	10000
seg_003	contig_2	1	3500
```

**ALN File (`ifn_aln`):**

Binary alignment file created by the `construct` command.

## Output File Formats

The command always writes `segment_summary.txt`, `adjacency_matrix.txt`, and `reach_matrix.txt`. With `-output_read_details T`, it also writes per-read association tables. With `-ifn_segment_clusters`, it also writes cluster-level matrix and summary files (see Parameters).

### segment_summary.txt

Summary statistics for each segment.

**Columns:**
* `segment`: Segment identifier
* `contig`: Contig name
* `start`: Segment start position (1-based inclusive)
* `end`: Segment end position (1-based inclusive)
* `length`: Segment length in bp
* `is_short`: `T` if segment length is less than `side_length + side_margin`; short segments use degenerate side coordinates and are skipped when creating associations
* `total_read_count_left`: Number of reads that enter/exit the segment on the left side
* `total_read_count_right`: Number of reads that enter/exit the segment on the right side
* `associated_read_count_left`: Reads that produced at least one association on the left side
* `associated_read_count_right`: Reads that produced at least one association on the right side
* `associated_segment_count_left`: Unique segments associated on the left side
* `associated_segment_count_right`: Unique segments associated on the right side

### adjacency_matrix.txt

Counts reads supporting direct adjacency between segment sides within the distance threshold.

**Columns:**
* `seg_src`, `seg_tgt`: Segment identifiers (ordered pair along the read)
* `side_src`, `side_tgt`: `"left"` or `"right"` for each segment
* `contig_src`, `start_src`, `end_src`: Segment coordinates for `seg_src` (1-based inclusive start/end)
* `contig_tgt`, `start_tgt`, `end_tgt`: Segment coordinates for `seg_tgt`
* `total_read_count`: Marginal total read count for `seg_src` on `side_src`
* `associated_read_count`: Marginal associated read count for `seg_src` on `side_src`
* `count`: Number of reads supporting this adjacency (read-space gap ≤ `max_adjacency_distance`)

For each `seg_src` ≠ `seg_tgt` entry, a second symmetric row is emitted with source and target swapped (same `count`, marginals taken from the swapped segment/side).

### reach_matrix.txt

Counts reads supporting connections between segment sides without adjacency constraint.

**Columns:** Same header and layout as `adjacency_matrix.txt`; `count` is the number of reads supporting the connection with no gap constraint.

## Matrix Definitions

* **Adjacency matrix**: counts ordered pairs of segment sides where the read exits one segment and enters the next with the estimated read-space gap between the intervals ≤ `max_adjacency_distance`. Each row tallies how many reads support that tight, locally consistent junction.
* **Reach matrix**: counts the same ordered pairs without the gap constraint, capturing longer-range connections along a read (potentially across multiple intermediate segments). Every adjacency is also a reach entry, but reach additionally records distant hops.

## Algorithm

1. **Segment side definition**: Let `min_segment_length = side_length + side_margin`. Segments with `length >= min_segment_length` define two side regions (0-based half-open coordinates on the contig):
   * **Left side:** the `side_length` bp interval starting `side_margin` bp after the segment start.
   * **Right side:** the `side_length` bp interval ending `side_margin` bp before the segment end (clamped so intervals stay inside the segment).
   Segments shorter than `min_segment_length` are marked short: side intervals collapse to the segment ends and are skipped when creating associations.

2. **Read processing**:
   * For each read, get all alignments sorted by read start
   * Apply parsimony filtering to select non-overlapping alignments:
     * Filter out alignments shorter than `side_length` (in contig coordinates)
     * Filter out alignments with mutation percentage > `max_mutation_percent` (over full alignment span; see parameter definition above)
     * Select parsimony subset: iteratively choose alignments prioritizing contigs with most selected read coverage, taking longest alignments (in read coordinates) that don't overlap previously selected alignments by more than `max_margin` (in read coordinates)
   * Intersect each alignment with segments to create intervals
   * For each interval, compute:
     * Side coverage (whether the read fully covers the shifted left/right side regions in contig coordinates, subject to the indel bp check below)
     * Mutation percentage within each side region (same event-count definition as above, excluding short indels)
   * Before marking a side as covered, check that the total indel bp within **that side’s interval only** (deletions overlapping it + insertions whose anchor position lies inside it, each with length ≥ `min_indel_length`) does not exceed `max_side_indel_bp`. This guards against large indels that overlap the side window itself. Indels entirely outside the left/right side intervals do not contribute to this sum (see **Limitations** below).

3. **Parsimony alignment selection**:
   * Quality filtering: remove alignments with contig span < `side_length` or mutation percentage > `max_mutation_percent`
   * Maintain selected alignments and track total read span per contig
   * Iteratively select alignments:
     * Choose top alignment from contigs with highest total selected read coverage (longest alignment in read coordinates per contig)
     * If no candidates remain for selected contigs, choose globally longest alignment (in read coordinates)
     * Accept alignment if read-coordinate overlap with all previously selected alignments ≤ `max_margin`
     * Update contig read coverage tracking
   * Tie-breaking for equally long alignments (in read coordinates):
     1. Prefer longer contig span
     2. Prefer fewer mutations
     3. Prefer lower read start coordinate
     4. Prefer lower contig index

4. **Interval sequence analysis**:
   * Traverse the interval sequence in read order
   * For every interval pointing in (left side covered on plus strand, or right side on minus strand) with mutation percentage ≤ `max_mutation_percent` on the relevant side, create reach/adjacency associations with all previously observed exit intervals
   * After processing entries, record new exit intervals (right side covered on plus strand, or left side on minus strand) whose mutation percentage is ≤ `max_mutation_percent`

5. **Adjacency check**:
   * Associations are considered adjacent when the estimated gap between their read coordinates is ≤ `max_adjacency_distance` bp.

## Examples

**Basic usage:**

```bash
alntools seg_matrix -ifn_aln output/sample.aln \
                    -ifn_segments segments.txt \
                    -odir output/seg_matrix
```

**Custom parameters:**

```bash
alntools seg_matrix -ifn_aln output/sample.aln \
                    -ifn_segments segments.txt \
                    -odir output/seg_matrix \
                    -max_mutation_percent 0.5 \
                    -max_adjacency_distance 500 \
                    -max_margin 20 \
                    -min_indel_length 3 \
                    -max_side_indel_bp 200 \
                    -side_length 800 \
                    -side_margin 150
```

## Notes

* **Short segments:** segments strictly shorter than `side_length + side_margin` are reported in the summary but skipped when creating reach or adjacency associations.
* **Side coverage:** the alignment must span the full left or right side interval in contig coordinates before that side is considered covered for entry/exit logic.
* **Mutation percentage on a side:** computed only over that side’s interval; mutations outside that interval do not affect that side’s percentage.
* **`max_side_indel_bp`:** uses the same `min_indel_length` threshold: only indels at or above that length contribute. For deletions, the contribution is the number of deleted bases overlapping the side window (a deletion starting before the window can still partially or fully consume it). For insertions, the contribution is the full insertion length when the insertion anchor lies inside the side window.

## Limitations

* **`max_side_indel_bp` is local to each side window.** Indels in the segment interior, in the margin strip before the left side window, after the right side window, or otherwise outside the defined left/right intervals are **not** added to the indel-bp total for that side. A read can still be counted as covering a side if it spans that interval in reference coordinates and indels elsewhere are arbitrarily large.
* **Global alignment filtering** uses the same event-count mutation percentage over the full alignment span; a single very large insertion is still one event, so it may not trigger `max_mutation_percent` unless the alignment is short or other mutations are dense.
