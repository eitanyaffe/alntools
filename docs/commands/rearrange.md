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
* `-max_margin <int>`: Maximum margin tolerance of alignments required to be adjacent (default: 10)
* `-min_element_length <int>`: Minimum element length for deletions and insertions (default: 50)
* `-min_anchor_length <int>`: Minimum anchor alignment length (default: 200)
* `-max_mutations_percent <double>`: Maximum mutations percentage for anchors (default: 0.1)

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
- Sequences that connect the alignments within the read, in case there is a gap between alignments or they overlap 
- **Read seams**: read sequences between adjacent alignments
- **Assembly seams**: assembly sequences between adjacent alignments

### From Contig to Read: The Mutation Process

Each rearrangement type follows a specific pattern for combining alignments and seams to reconstruct the read sequence:

#### Large Insertion
**Pattern:** A → left_seam → X → right_seam → B

1. **Clip at boundaries**: Remove sequences outside the rearrangement region
2. **Handle gaps/overlaps**: If there's a gap between A and the insertion point, keep it in the assembly seam; overlaps are also preserved in seams
3. **Add left seam**: Bridge from anchor A to element X
4. **Add element X**: Insert the element sequence (properly oriented)  
5. **Add right seam**: Bridge from element X to anchor B
6. **Mutations applied**: Each alignment (A, X, B) has its mutations applied during sequence extraction

#### Large Inversion  
**Pattern:** A → left_seam → X_reversed → right_seam → B

Same process as insertion, but the element X is reverse-complemented. The element sequence is right there in the middle, just flipped in orientation relative to the flanking anchors.

#### Large Deletion
**Pattern:** A → seam → B (skipping deleted region)

1. **Clip at boundaries**: Remove sequences outside the rearrangement region
2. **Skip deleted region**: Jump directly from anchor A to anchor B
3. **Add bridging seam**: Connect A and B with any intervening read sequence
4. **No element**: The deleted assembly region is absent from the final read sequence

**Note:** These instructions deal only with rearrangement of alignments and seams. To get the actual read sequence, alignments are also mutated along the way as needed using their individual mutation profiles.

### Seams Explained

**Seams represent the boundaries between alignments:**
- **Gap seams**: When there's a gap in the read, the seam contains the bridging sequence
- **Overlap seams**: When alignments overlap, the seam represents the overlapping region
- Seams are essential for reconstructing the complete read sequence from individual alignments

**Seam Format in Output:**
- **Gap seams**: Prefixed with `+` (e.g., `+ATCG` indicates a 4bp gap with sequence ATCG)
- **Overlap seams**: Prefixed with `-` (e.g., `-GCTA` indicates a 4bp overlap with sequence GCTA)
- **Multiple seams**: Separated by colons (`:`) for events with multiple seam regions
- **Empty seams**: No prefix, represented as empty strings between colons

**Examples:**
- `+ATCG` - Single gap seam with sequence ATCG
- `-GC` - Single overlap seam with sequence GC
- `+ATCG:+TTAA` - Two gap seams (for insertions/inversions with left and right seams)
- `+ATCG:` - Gap seam followed by empty seam
- `:+TTAA` - Empty seam followed by gap seam

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
