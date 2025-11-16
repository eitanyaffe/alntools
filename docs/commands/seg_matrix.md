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
* `-max_mutation_percent <double>`: Maximum mutation percentage allowed for side coverage (default: 1.0). Sides with higher mutation rates are ignored.
* `-max_adjacency_distance <int>`: Maximum distance in bp for adjacency (default: 200). Distance is measured in estimated read space.
* `-max_margin <int>`: Maximum overlap between consecutive alignments in bp (default: 10). Alignments overlapping more than this threshold are skipped.
* `-min_indel_length <int>`: Minimum indel length to include in mutation density calculations (default: 3). Shorter indels are excluded.
* `-max_mutations_percent <double>`: Maximum mutation percentage for alignment filtering (default: 1.0). Alignments exceeding this threshold are skipped entirely.
* `-side_length <int>`: Length in bp of each segment side used for coverage analysis (default: 1000).
* `-side_margin <int>`: Margin in bp between the segment edge and the start of each side region (default: 200).

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

The command generates three tab-delimited files in the output directory:

### segment_summary.txt

Summary statistics for each segment.

**Columns:**
* `segment`: Segment identifier
* `contig`: Contig name
* `start`: Segment start position (1-based inclusive)
* `end`: Segment end position (1-based inclusive)
* `length`: Segment length in bp
* `is_short`: Flag indicating if segment is shorter than `2 * (side_length + side_margin)` (T if short, F otherwise)
* `total_read_count_left`: Number of reads that enter/exit the segment on the left side
* `total_read_count_right`: Number of reads that enter/exit the segment on the right side
* `associated_read_count_left`: Reads that produced at least one association on the left side
* `associated_read_count_right`: Reads that produced at least one association on the right side
* `associated_segment_count_left`: Unique segments associated on the left side
* `associated_segment_count_right`: Unique segments associated on the right side

### adjacency_matrix.txt

Counts reads supporting direct adjacency between segment sides within the distance threshold.

**Columns:**
* `seg1`: First segment identifier
* `seg2`: Second segment identifier
* `side1`: Side of seg1 ("left" or "right")
* `side2`: Side of seg2 ("left" or "right")
* `contig1`, `start1`, `end1`: Coordinates for seg1
* `contig2`, `start2`, `end2`: Coordinates for seg2
* `total_read_count_1`, `associated_read_count_1`: Marginal counts for seg1 on `side1`
* `total_read_count_2`, `associated_read_count_2`: Marginal counts for seg2 on `side2`
* `count`: Number of reads supporting this adjacency (read-space gap ≤ `max_adjacency_distance`)

### reach_matrix.txt

Counts reads supporting connections between segment sides without adjacency constraint.

**Columns:**
* `seg1`: First segment identifier
* `seg2`: Second segment identifier
* `side1`: Side of seg1 ("left" or "right")
* `side2`: Side of seg2 ("left" or "right")
* `contig1`, `start1`, `end1`: Coordinates for seg1
* `contig2`, `start2`, `end2`: Coordinates for seg2
* `total_read_count_1`, `associated_read_count_1`: Marginal counts for seg1 on `side1`
* `total_read_count_2`, `associated_read_count_2`: Marginal counts for seg2 on `side2`
* `count`: Number of reads supporting this connection

## Matrix Definitions

* **Adjacency matrix**: counts ordered pairs of segment sides where the read exits one segment and enters the next with the estimated read-space gap between the intervals ≤ `max_adjacency_distance`. Each row tallies how many reads support that tight, locally consistent junction.
* **Reach matrix**: counts the same ordered pairs without the gap constraint, capturing longer-range connections along a read (potentially across multiple intermediate segments). Every adjacency is also a reach entry, but reach additionally records distant hops.

## Algorithm

1. **Segment Side Definition**: Each segment that is at least `2 * (side_length + side_margin)` bp long defines two side regions:
   * Left side: the `side_length` bp interval starting `side_margin` bp from the segment start
   * Right side: the `side_length` bp interval ending `side_margin` bp before the segment end
   Segments shorter than this threshold are marked as short segments and skipped when creating associations.

2. **Read Processing**:
   * For each read, get all alignments sorted by read start
   * Filter alignments with overlap > `max_margin`
   * Intersect each alignment with segments to create intervals
   * For each interval, compute:
     * Side coverage (whether the read fully covers the shifted left/right side regions)
     * Mutations within each covered side (excluding short indels)

3. **Interval Sequence Analysis**:
   * Traverse the interval sequence in read order
   * For every interval pointing in (left side covered on plus strand, or right side on minus strand) with mutation percentage ≤ `max_mutation_percent`, create reach/adjacency associations with all previously observed exit intervals
   * After processing entries, record new exit intervals (right side covered on plus strand, or left side on minus strand) whose mutation percentage is ≤ `max_mutation_percent`

4. **Adjacency Check**:
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
                    -side_length 800 \
                    -side_margin 150
```

## Notes

* Segments shorter than `2 * (side_length + side_margin)` bp are marked as short segments. They are reported in the summary table but skipped when creating reach or adjacency associations.
* Side coverage requires the read to fully span the shifted side regions (including the side_margin offset) before a segment is considered to be pointing in or out.
* Mutation percentages used for side coverage are computed ignoring indels shorter than `min_indel_length`.