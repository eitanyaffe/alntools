# rearrange

Detects genome rearrangements (large insertions, deletions, and inversions) from read alignments by analyzing alignment patterns within reads.

## Syntax

```bash
alntools rearrange -ifn_aln <input.aln> -ofn_prefix <output_prefix> [options]
# OR for multi-library analysis:
alntools rearrange -ifn_libraries <libraries.tsv> -ofn_prefix <output_prefix> [options]
```

## Parameters

**Mandatory Arguments:**
* Either `-ifn_aln <fn>`: Single ALN file (treated as library "sample")
* Or `-ifn_libraries <fn>`: Tab-delimited file with library definitions (format: `id fn`)
* `-ofn_prefix <fn>`: Output prefix for result files

**Optional Arguments:**
* `-ifn_intervals <fn>`: Tab-delimited file with query intervals (format: `contig start end`). Only process rearrangements where anchor alignments overlap these intervals
* `-max_margin <int>`: Maximum margin tolerance between anchor alignments (default: 10)
* `-min_element_length <int>`: Minimum element length for deletions and insertions (default: 50)
* `-min_anchor_length <int>`: Minimum anchor alignment length (default: 200)
* `-max_mutations_percent <double>`: Maximum mutations percentage for all alignments (default: 0.01)

## Examples

```bash
# Single library rearrangement detection
alntools rearrange -ifn_aln output/sample.aln -ofn_prefix output/rearrangements

# Multi-library analysis with interval filtering
echo -e "id\tfn\nlib1\toutput/sample1.aln\nlib2\toutput/sample2.aln" > libraries.tsv
echo -e "contig\tstart\tend\nchr1\t1000000\t2000000" > intervals.tsv
alntools rearrange -ifn_libraries libraries.tsv -ifn_intervals intervals.tsv -ofn_prefix output/rearrangements
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command
- **Libraries file**: Tab-delimited file with columns `id` and `fn`
- **Intervals file**: See [File Formats](../file-formats.md#intervals-format)

## Rearrangement Types

- **large_insert**: Large insertions where element sequence comes from another location
- **large_delete**: Large deletions with gaps between anchor alignments  
- **large_invert**: Large inversions where element sequence is in reverse orientation

## How Rearrangement Events Work

Rearrangement detection analyzes reads with multiple alignments to identify structural variations. Each rearrangement event is defined by:

### Event Components

**Anchor Alignments (A and B):**
- Two high-quality alignments that flank the rearrangement
- Alignment A: upstream anchor
- Alignment B: downstream anchor
- Both must meet minimum length and mutation rate criteria

**Element Alignment (X, optional):**
- For insertions and inversions: alignment representing the inserted/inverted sequence
- For deletions: no element alignment (gap between anchors)

**Seams:**
- Sequences that connect the alignments within the read
- **Read seams**: actual sequences from the read that bridge alignments
- **Assembly seams**: corresponding sequences from the assembly (for comparison)

### From Contig to Read: The Mutation Process

To understand how a rearrangement transforms a contig segment into a read segment:

1. **Start with anchor A**: Extract contig sequence and apply mutations from alignment A
2. **Add left seam**: Insert read sequence (for insertions) or skip sequence (for deletions)
3. **Add element X** (if present): Extract and mutate element sequence, applying reverse complement for inversions
4. **Add right seam**: Insert read sequence (for insertions) or skip sequence (for deletions)  
5. **End with anchor B**: Extract contig sequence and apply mutations from alignment B

**For deletions:** A → left_seam → B (skipping deleted region)
**For insertions:** A → left_seam → X → right_seam → B (inserting element X)
**For inversions:** A → left_seam → X_reversed → right_seam → B (reversing element X)

### Seams Explained

**Seams represent the boundaries between alignments:**
- **Gap seams**: When there's a gap in the read, the seam contains the bridging sequence
- **Overlap seams**: When alignments overlap, the seam represents the overlapping region
- Seams are essential for reconstructing the complete read sequence from individual alignments

**To get seams, you need:**
- Read sequences (to extract actual bridging sequences)
- Alignment coordinates (to identify gaps and overlaps)
- Assembly sequences (to compare expected vs. observed seams)

## Verification Process

The verification component checks consistency between detected events and the original sequences:

**Purpose:** Ensure that the rearrangement model accurately explains how the contig segments were mutated to produce the observed read sequence.

**Process:**
1. **Extract sequences**: Get contig and read sequences for all alignments
2. **Apply mutations**: Mutate contig segments according to alignment mutations
3. **Construct final sequence**: Assemble A + seams + X + seams + B
4. **Compare**: Verify that the constructed sequence matches the actual read sequence

**Verification is not about finding rearrangements** - it's about validating that the detected rearrangements are consistent with the underlying sequence data.

## Output File Formats

**Single Library Mode:**
- `{prefix}_sample_read_events.tsv`: Individual rearrangement events per read
- `{prefix}_sample_aggregated_events.tsv`: Aggregated events with read counts and coverage
- `{prefix}_sample_read_support.tsv`: Reads supporting each event (including coverage reads)

**Multi-Library Mode (additional files):**
- `{prefix}_events.tsv`: Combined events across all libraries with library counts
- `{prefix}_support.tsv`: Read support matrix (events × libraries)  
- `{prefix}_coverage.tsv`: Coverage matrix (events × libraries)

See [File Formats](../file-formats.md) for detailed specifications of output formats.

## Notes

- Rearrangement detection requires reads with multiple alignments to the same contig
- The algorithm is designed to detect large-scale structural variations (>50bp by default)
- Verification ensures that detected events are biologically plausible
- Multi-library analysis enables population-level structural variation studies
- Use interval filtering to focus on specific genomic regions of interest
