# File Formats

This document details the file formats used by `alntools` for both input and output.

Positions are 1-based in input and output files, and represented internally using a 0-based coordinate system.

## Input Formats

### PAF Format

`alntools` processes alignment data in the Pairwise mApping Format (PAF). The following columns are required:

| Column Index | Description                                      | Type  |
|--------------|--------------------------------------------------|-------|
| 1            | Query sequence name                              | string|
| 2            | Query sequence length                            | int   |
| 3            | Query start                                      | int   |
| 4            | Query end                                        | int   |
| 5            | Relative strand ('+' or '-')                     | char  |
| 6            | Target sequence name (contig ID)                 | string|
| 7            | Target sequence length                           | int   |
| 8            | Target start on original strand                  | int   |
| 9            | Target end on original strand                    | int   |

Additionally, the **`cs:Z:` tag** (difference string) **must be present** in one of the optional fields (column 13 onwards). This tag encodes the base differences between the query (read) and the target (contig) and is essential for mutation analysis.

**Example PAF Line:**
```text
read_298    150 10  145 +   contig_12   5000    1500    1635    120 135 60  cs:Z::10*at:5+gg:20-c:15*cg:70
```

### Intervals Format

The intervals file is a tab-delimited file specifying regions to query:

| Column  | Description                           | Type  |
|---------|---------------------------------------|-------|
| contig  | Name of the contig                    | string|
| start   | Start position                        | int   |
| end     | End position                          | int   |

**Example intervals file:**
```
contig	start	end
ctg25860	120448	130902
ctg26175	230065	252386
```

### Libraries Format

For multi-library analysis (variants and rearrange commands):

| Column | Description                    | Type  |
|--------|--------------------------------|-------|
| id     | Library identifier             | string|
| fn     | Path to ALN file               | string|

**Example libraries file:**
```
id	fn
lib1	output/sample1.aln
lib2	output/sample2.aln
lib3	output/sample3.aln
```

### Gene Table Format

For gene annotation (optional):

| Column | Description                    | Type  |
|--------|--------------------------------|-------|
| gene   | Gene identifier                | string|
| contig | Contig name                    | string|
| start  | Gene start position            | int   |
| end    | Gene end position              | int   |
| strand | Strand (+ or -)                | char  |
| desc   | Gene description (optional)    | string|

## Output Formats

### Query Mode Outputs

#### Full Mode Output

**`{prefix}_alignments.tsv`:**

| Column         | Description                              | Type    |
|----------------|------------------------------------------|---------|
| alignment_index| Unique index for the alignment           | int     |
| read_id        | ID of the read                           | string  |
| contig_id      | ID of the contig                         | string  |
| read_start     | Start position on read                   | int     |
| read_end       | End position on read                     | int     |
| contig_start   | Start position on contig                 | int     |
| contig_end     | End position on contig                   | int     |
| is_reverse     | Whether alignment is on reverse strand   | boolean |
| cs_tag         | CIGAR string encoding differences        | string  |
| num_mutations  | Number of mutations in the alignment     | int     |
| height         | Vertical position for visualization      | int     |
| chunk_id       | Chunk identifier for height calculation  | int     |

**`{prefix}_mutations.tsv`:**

| Column         | Description                              | Type    |
|----------------|------------------------------------------|---------|
| alignment_index| Index of parent alignment                | int     |
| read_id        | ID of the read                           | string  |
| contig_id      | ID of the contig                         | string  |
| type           | Mutation type (SUB/INS/DEL)              | string  |
| position       | Position on contig                       | int     |
| desc           | Description of the mutation              | string  |
| height         | Vertical position for visualization      | int     |

**`{prefix}_reads.tsv`:**

| Column               | Description                              | Type    |
|----------------------|------------------------------------------|---------|
| read_id              | ID of the read                           | string  |
| contig_id            | ID of the contig                         | string  |
| read_length          | Total length of the read                 | int     |
| span_start           | Start of read's aligned region on contig| int     |
| span_end             | End of read's aligned region on contig  | int     |
| total_aligned_length | Total aligned base pairs                 | int     |
| num_alignments       | Number of alignments for this read      | int     |
| num_mutations        | Total mutations across all alignments   | int     |
| height               | Vertical position for visualization      | int     |

**`{prefix}_chunks.tsv`:**

| Column               | Description                              | Type    |
|----------------------|------------------------------------------|---------|
| chunk_id             | Unique chunk identifier                  | int     |
| read_id              | ID of the read                           | string  |
| contig_id            | ID of the contig                         | string  |
| chunk_start          | Start position of chunk on contig       | int     |
| chunk_end            | End position of chunk on contig         | int     |
| num_alignments       | Number of alignments in chunk            | int     |
| num_mutations        | Total mutations in chunk                 | int     |
| height               | Vertical position for visualization      | int     |

#### Pileup Mode Output

**`{prefix}_pileup.tsv`:**

| Column   | Description                                 | Type   |
|----------|---------------------------------------------|--------|
| contig   | Contig ID                                   | string |
| position | Position on contig                          | int    |
| variant  | Variant type (REF or mutation description)  | string |
| count    | Count of occurrences                        | int    |
| coverage | Total read coverage at this position        | int    |
| cumsum   | Cumulative count                            | int    |

#### Bin Mode Output

**`{prefix}_bins.tsv`:**

| Column                | Description                               | Type   |
|-----------------------|-------------------------------------------|--------|
| contig                | Contig ID                                 | string |
| bin_start             | Bin start position                        | int    |
| bin_end               | Bin end position                          | int    |
| bin_length            | Bin length                                | int    |
| sequenced_bp          | Sequenced base pairs in bin               | int    |
| mutation_count        | Number of mutations in bin                | int    |
| seg_sites_density     | Segregating sites density (sites per bp) | double |
| non_ref_sites_density | Non-reference sites density (sites per bp)| double |
| dist_none             | Fraction with 0 mutations                 | double |
| dist_5                | Fraction with 1e-5 to 1e-4 mutations/bp  | double |
| dist_4                | Fraction with 1e-4 to 1e-3 mutations/bp  | double |
| dist_3                | Fraction with 1e-3 to 1e-2 mutations/bp  | double |
| dist_2                | Fraction with 1e-2 to 1e-1 mutations/bp  | double |
| dist_1_plus           | Fraction with >1e-1 mutations/bp         | double |

#### Consensus Mode Output

**`{prefix}_consensus.tsv`:**

| Column       | Description                              | Type   |
|--------------|------------------------------------------|--------|
| contig       | Contig ID                                | string |
| position     | Position on contig                       | int    |
| variant_type | Variant type (SUB/INS/DEL)               | string |
| variant_desc | Variant description (A:G, +ATG, -CC)     | string |
| count        | Supporting reads                         | int    |
| coverage     | Total reads at position                  | int    |
| frequency    | Variant frequency (count/coverage)       | double |

#### Variants Mode Output

**`{prefix}_variants.tsv`:**

| Column         | Description                              | Type   |
|----------------|------------------------------------------|--------|
| variant_id     | Sequential variant ID (v1, v2, v3...)   | string |
| contig         | Contig name                              | string |
| coord          | 1-based coordinate                       | int    |
| type           | Variant type (sub/ins/del/left_clip/right_clip) | string |
| sequence       | Variant sequence                         | string |
| desc           | Human-readable description               | string |
| library_count  | Number of libraries containing variant   | int    |
| total_support  | Total supporting reads across libraries  | int    |
| total_coverage | Total coverage across libraries          | int    |
| frequency      | Total support / total coverage ratio     | double |

**`{prefix}_support.tsv`:**
Matrix with variant IDs as rows and library IDs as columns (read support counts)

**`{prefix}_coverage.tsv`:**
Matrix with variant IDs as rows and library IDs as columns (coverage counts)

### Coverage Mode Output

**`{prefix}_coverage.txt`:**

| Column                | Description                              | Type |
|-----------------------|------------------------------------------|------|
| contig_count          | Number of contigs                        | int  |
| read_count            | Number of reads                          | int  |
| alignment_count       | Number of alignments                     | int  |
| total_read_bp         | Total read base pairs                    | int  |
| total_assembly_bp     | Total assembly base pairs                | int  |
| contig_n50            | N50 statistic of contig lengths          | int  |
| aligned_count         | Number of reads with alignments          | int  |
| non_aligned_read_bp   | Base pairs in unaligned read regions    | int  |
| non_aligned_contig_bp | Base pairs in unaligned contig regions  | int  |

**`{prefix}_non_aligned_reads.txt`:**

| Column  | Description                    | Type   |
|---------|--------------------------------|--------|
| read_id | Read identifier                | string |
| start   | Start position (1-based)       | int    |
| end     | End position (1-based)         | int    |
| length  | Length of non-aligned region   | int    |

**`{prefix}_non_aligned_contigs.txt`:**

| Column    | Description                    | Type   |
|-----------|--------------------------------|--------|
| contig_id | Contig identifier              | string |
| start     | Start position (1-based)       | int    |
| end       | End position (1-based)         | int    |
| length    | Length of non-aligned region   | int    |

### Breaks Mode Output

**`{output}.tsv`:**

| Column      | Description                              | Type   |
|-------------|------------------------------------------|--------|
| contig      | Contig identifier                        | string |
| position    | 1-based position on contig               | int    |
| orientation | "left" (read starts) or "right" (ends)  | string |
| t           | Observed reads at position               | int    |
| e           | Expected reads per position              | double |
| enrichment  | Observed/expected ratio (t/e)            | double |
| pval        | Raw p-value from binomial test           | double |
| qval        | Benjamini-Hochberg corrected q-value     | double |

### Rearrange Mode Output

**Single Library Mode:**

**`{prefix}_sample_read_events.tsv`:**

| Column           | Description                              | Type   |
|------------------|------------------------------------------|--------|
| read_id          | Read identifier                          | string |
| event_id         | Event identifier                         | string |
| type             | Rearrangement type                       | string |
| contig_id        | Contig of anchor alignments              | string |
| read_strand      | Read orientation (+ or -)                | string |
| out_clip         | Outer clip position (contig coords)      | int    |
| in_clip          | Inner clip position (contig coords)      | int    |
| read_clip_out    | Outer clip position (read coords)        | int    |
| read_clip_in     | Inner clip position (read coords)        | int    |
| span_start       | Start of event span (contig coords)      | int    |
| span_end         | End of event span (contig coords)        | int    |
| read_span_start  | Start of event span (read coords)        | int    |
| read_span_end    | End of event span (read coords)          | int    |
| element_contig   | Element contig (empty for deletions)     | string |
| element_strand   | Element strand (empty for deletions)     | string |
| element_start    | Element start (0 for deletions)          | int    |
| element_end      | Element end (0 for deletions)            | int    |
| read_seams       | Read seam sequences (see format below)   | string |
| assembly_seams   | Assembly seam sequences (see format below) | string |

**Seam Sequence Format:**

Seam sequences in both `read_seams` and `assembly_seams` columns use the following format:
- **Multiple seams**: Separated by colons (`:`)
- **Gap seams**: Prefixed with `+` (e.g., `+ATCG` indicates a 4bp gap with sequence ATCG)
- **Overlap seams**: Prefixed with `-` (e.g., `-GCTA` indicates a 4bp overlap with sequence GCTA)
- **Empty seams**: No prefix, empty string between colons

**`{prefix}_sample_aggregated_events.tsv`:**
Aggregated version with read counts and coverage statistics per event.

**`{prefix}_sample_read_support.tsv`:**
Reads supporting each aggregated event.

**Multi-Library Mode (additional files):**

**`{prefix}_events.tsv`:**
Combined events across libraries with library counts.

**`{prefix}_support.tsv`:**
Read support matrix (events × libraries).

**`{prefix}_coverage.tsv`:**
Coverage matrix (events × libraries).

### Gene Annotation Output

When gene annotation is enabled, additional columns and files are generated:

**Additional columns in main output files:**
- `is_genic`: Boolean indicating if variant is within a gene

**`{prefix}_genic.tsv`:**

| Column        | Description                              | Type   |
|---------------|------------------------------------------|--------|
| row_id        | Row identifier from main output          | int    |
| gene_id       | Gene identifier                          | string |
| gene_desc     | Gene description                         | string |
| aa_coord      | Amino acid coordinate                    | int    |
| variant_codon | Variant codon                            | string |
| ref_codon     | Reference codon                          | string |
| variant_type  | Mutation type (syn/non-syn)              | string |
| mutation_desc | Mutation description (S:F, syn, etc.)    | string |

**`{prefix}_intergenic.tsv`:**

| Column           | Description                              | Type   |
|------------------|------------------------------------------|--------|
| row_id           | Row identifier from main output          | int    |
| gene_left        | Left flanking gene                       | string |
| gene_right       | Right flanking gene                      | string |
| orientation_left | Left gene orientation                    | string |
| orientation_right| Right gene orientation                   | string |
| distance_left    | Distance to left gene                    | int    |
| distance_right   | Distance to right gene                   | int    |

## Binary Formats

### ALN Format

The ALN format is a binary format optimized for fast querying. It stores:
- Read and contig identifiers
- Alignment coordinates and orientations
- Mutation information from PAF cs:Z: tags
- Indexing structures for efficient queries

Use the `extract` command to convert ALN files to readable tab-delimited formats for inspection.
