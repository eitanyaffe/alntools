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
- **Assembly seams**: assembly sequences that define how alignments interact in assembly space

### From Contig to Read: The Mutation Process

Each rearrangement type follows a specific pattern for combining mutated alignment segments and seams to reconstruct the read sequence:

#### Large Insertion
**Pattern:** A → left_seam → X → right_seam → B

1. **Create A sequence**: Extract contig sequence from alignment A coordinates, apply A's mutations
2. **Handle left seam**: 
   - If gap seam: add the bridging sequence from the read
   - If overlap seam: verify exact match and remove overlapping sequence from end of A sequence
3. **Create X sequence**: Extract element sequence from alignment X coordinates, apply X's mutations
4. **Apply orientation**: If X alignment is reverse strand, reverse-complement the X sequence
5. **Handle right seam**: Same gap/overlap logic as left seam
6. **Create B sequence**: Extract contig sequence from alignment B coordinates, apply B's mutations

#### Large Inversion  
**Pattern:** A → left_seam → X_reversed → right_seam → B

Same process as insertion. The element X sequence is reverse-complemented if the X alignment is on the reverse strand, creating the inversion effect.

#### Large Deletion
**Pattern:** A → seam → B

1. **Create A sequence**: Extract and mutate alignment A sequence
2. **Handle seam**: 
   - If gap seam: add bridging sequence from read
   - If overlap seam: verify exact match and remove overlapping sequence from A sequence
3. **Create B sequence**: Extract and mutate alignment B sequence
4. **No element**: The deleted assembly region between A and B is skipped entirely

**Key Points:**
- **Clip coordinates** (`out_clip`, `in_clip`) define the breakpoints: `out_clip = min(A.end, B.end)`, `in_clip = max(A.start, B.start)`
- **Mutations are applied** to each alignment segment (A, B, X) using their individual PAF mutation profiles
- **Seam handling** depends on whether alignments have gaps (add sequence) or overlaps (remove sequence) between them

### Assembly Seams: How Alignments Interact in Assembly Space

Assembly seams define the relationship between alignments in the reference assembly coordinates. They are not used for read construction but describe the structural variation in assembly space:

#### Large Insertion
**Assembly pattern:** A ← seam → B (one seam between A and B)

The assembly seam represents the sequence between anchor alignments A and B in the assembly. This shows what assembly sequence is "replaced" by the inserted element X in the read.

#### Large Deletion  
**Assembly pattern:** A → B (no seam)

No assembly seam because the deleted region is the gap between A and B alignments in the assembly. The deletion is simply the missing assembly sequence between the anchors.

#### Large Inversion
**Assembly pattern:** A ← left_seam → X ← right_seam → B (two seams)

- **Left seam (A-X)**: Assembly sequence between anchor A and element X
- **Right seam (X-B)**: Assembly sequence between element X and anchor B

This shows how the element X relates to the flanking anchors A and B in the original assembly coordinates.

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
