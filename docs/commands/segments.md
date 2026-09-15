# segments

Performs segmentation analysis by detecting breakpoints in read alignments and generating genomic segments based on dangling read patterns across libraries.

## Syntax

```bash
alntools segments -ifn_aln <input.aln> -odir <output_dir> [options]
# OR for multi-library analysis:
alntools segments -ifn_libraries <libraries.txt> -odir <output_dir> [options]
```

## Parameters

**Mandatory Arguments:**
* Either `-ifn_aln <fn>`: Single ALN file (treated as library "sample")
* Or `-ifn_libraries <fn>`: Tab-delimited file with library definitions (format: `id fn`)
* `-odir <dir>`: Output directory for result files

**Breakpoint Detection Parameters:**
* `-max_margin <int>`: Maximum margin tolerance for breakpoint clustering (default: 20)
* `-min_anchor_length <int>`: Minimum anchor alignment length (default: 1000)
* `-min_dangle_length <int>`: Minimum length for an alignment to count as a continuation of the anchor (default: 1000)
* `-max_anchor_mutations_percent <double>`: Maximum mutations percentage for anchor alignments (default: 0.1)
* `-min_alignment_distance <int>`: Minimum distance between anchor and a continuation on the same contig (default: 200)

**Breakpoint Selection Parameters:**
* `-min_breakpoint_read_support <int>`: Minimum read support for selecting breakpoints (default: 2)
* `-min_breakpoint_frequency <double>`: Minimum frequency for selecting breakpoints (default: 0.2)

**Segment Generation Parameters:**
* `-min_segment_length <int>`: Minimum segment length for filtering breakpoints (default: 200)
* `-max_segment_length <int>`: Segments longer than this are split by artificial breakpoints; 0 disables splitting (default: 0)

## Examples

```bash
# Single library segmentation analysis
alntools segments -ifn_aln output/sample.aln -odir output/segments

# Multi-library analysis with custom parameters
echo -e "id\tfn\nlib1\toutput/sample1.aln\nlib2\toutput/sample2.aln" > libraries.txt
alntools segments -ifn_libraries libraries.txt -odir output/segments \
  -min_segment_length 500 \
  -min_breakpoint_read_support 3 \
  -min_breakpoint_frequency 0.3
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command
- **Libraries file**: Tab-delimited file with columns `id` and `fn`

## How Segmentation Works

### 1. Breakpoint Detection

For each library, the tool identifies potential breakpoints by analyzing read alignment patterns:

A breakpoint is called at an anchor boundary when the read does **not** continue
along the same contig past that boundary, i.e. the read dangles away from the contig.

**Breakpoint Types:**
- **dangle_left**: Read dangles off the left side of the anchor (breakpoint at the anchor's contig start)
- **dangle_right**: Read dangles off the right side of the anchor (breakpoint at the anchor's contig end)

Strand is taken into account: for a reverse-strand anchor the read-order sense of the two
boundaries is swapped.

**Detection Criteria:**

Each alignment of the read is tested as an anchor, and each of its two boundaries is
tested separately:
- Anchor must be at least `min_anchor_length` bp with ≤ `max_anchor_mutations_percent` mutations
- The read must extend at least `min_alignment_distance + min_dangle_length` bp past the
  anchor boundary, so that a continuation could have been observed had one existed
- No **continuation** may exist on that side. An alignment counts as a continuation when it
  lies within `min_alignment_distance` bp of the anchor in read coordinates, is on the same
  contig and strand, is at least `min_dangle_length` bp long, and is within
  `min_alignment_distance` bp of the anchor in contig coordinates

Note that `min_dangle_length` and `min_alignment_distance` act purely as rejection
thresholds on candidate continuations; the dangling part of the read itself is not measured.

### 2. Breakpoint Aggregation

Breakpoints from all libraries are clustered into aggregate breakpoints:
- Breakpoints within `max_margin` bp on the same contig and type are grouped
- Coordinate is determined by the mode (most frequent position) within the cluster
- Support and coverage are calculated per library

### 3. Breakpoint Filtering

Aggregate breakpoints are selected based on two criteria (either passes):
- **Global threshold**: `read_support` ≥ `min_breakpoint_read_support` AND `frequency` ≥ `min_breakpoint_frequency`
- **Per-library threshold**: At least one library meets both thresholds independently

### 4. Segment Generation

Selected breakpoints are used to define genomic segments:
- Breakpoints are sorted by coordinate within each contig
- Only breakpoints that maintain segments ≥ `min_segment_length` are used
- Breakpoints too close to contig boundaries (< `min_segment_length`) are filtered out
- If `max_segment_length` > 0, any remaining gap longer than it is split evenly by
  artificial coordinates, which have no breakpoint id; splitting is disabled when it is 0
- Segments are created between consecutive accepted breakpoints

**Output:**
- All segments are retained (even if shorter than threshold)
- Warning is issued for segments shorter than `min_segment_length`

## Output Files

All output files are written to the specified output directory:

### read_breakpoints.txt
Individual breakpoint calls from reads, per library.

**Columns:**
- `lib_id`: Library identifier
- `read_id`: Read identifier
- `contig`: Contig name
- `coord`: Breakpoint coordinate (1-based)
- `type`: Breakpoint type (dangle_left/dangle_right)
- `anchor_length`: Length of anchor alignment
- `anchor_mutations`: Number of mutations in anchor
- `aggregate_breakpoint_id`: Id of the aggregate breakpoint this call was clustered into

### breakpoints.txt
Aggregated breakpoints across all libraries.

**Columns:**
- `breakpoint_id`: Unique breakpoint identifier (b1, b2, ...)
- `contig`: Contig name
- `coord`: Breakpoint coordinate (1-based)
- `type`: Breakpoint type (dangle_left/dangle_right)
- `read_support`: Total read support across all libraries
- `frequency`: Frequency (read_support / total_coverage)
- `selected`: Whether breakpoint passed selection thresholds (T/F)
- `is_segment_break`: Whether breakpoint was used to define segments (T/F)

### breakpoints_support.txt
Matrix of read support per breakpoint per library.

**Format:**
- Rows: Breakpoints
- Columns: Libraries
- Values: Number of supporting reads

### breakpoints_coverage.txt
Matrix of coverage per breakpoint per library.

**Format:**
- Rows: Breakpoints
- Columns: Libraries
- Values: Number of reads covering the breakpoint position

### segments.txt
Final genomic segments defined by filtered breakpoints.

**Columns:**
- `segment`: Segment identifier (s1, s2, ...)
- `contig`: Contig name
- `start`: Segment start coordinate (1-based, inclusive)
- `end`: Segment end coordinate (1-based, inclusive)
- `length`: Segment length in bp

## Notes

- Breakpoints marked as `selected=T` passed support/frequency thresholds
- Breakpoints marked as `is_segment_break=T` were actually used to create segments
- Some selected breakpoints may not be segment breaks if they violate `min_segment_length` constraints
- Segments shorter than `min_segment_length` may exist at contig boundaries; these generate a warning

