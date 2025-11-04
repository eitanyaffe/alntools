# cov_matrix

Computes per-base coverage and variance statistics for genomic segments across multiple libraries, producing a coverage matrix suitable for downstream binning analysis (e.g., MetaBAT).

## Syntax

```bash
alntools cov_matrix -ifn_libraries <library_table> \
                    -ifn_segments <segment_table> \
                    -ifn_fasta <contigs.fasta> \
                    -ofn_mat <coverage_matrix.txt> \
                    -ofn_fasta <segments.fasta>
```

## Parameters

**Mandatory Arguments:**
* `-ifn_libraries <fn>`: Library table with `lib_id` and `aln_fn` columns.
* `-ifn_segments <fn>`: Segment table with `segment`, `contig`, `start`, `end` columns.
* `-ifn_fasta <fn>`: Contig fasta file.
* `-ofn_mat <fn>`: Output coverage matrix file.
* `-ofn_fasta <fn>`: Output segment fasta file.

## Example

```bash
alntools cov_matrix -ifn_libraries libraries.txt \
                    -ifn_segments segments.txt \
                    -ifn_fasta contigs.fa \
                    -ofn_mat coverage_matrix.txt \
                    -ofn_fasta segments.fa
```

## Input File Formats

**Library Table (`ifn_libraries`):**

Tab-delimited file with header:
- `lib_id`: Library identifier
- `aln_fn`: Path to ALN file for this library

Example:
```
lib_id	aln_fn
early	output/early.aln
pre	output/pre.aln
post	output/post.aln
late	output/late.aln
```

**Segment Table (`ifn_segments`):**

Tab-delimited file with header (coordinates are 1-based, inclusive):
- `segment`: Unique segment identifier
- `contig`: Contig identifier
- `start`: Start position (1-based)
- `end`: End position (1-based, inclusive)

Example:
```
segment	contig	start	end
seg_001	contig_1	1	5000
seg_002	contig_1	5001	10000
seg_003	contig_2	1	3500
```

**Contig FASTA (`ifn_fasta`):**

Standard FASTA format containing the contig sequences referenced in the segment table.

## Output File Formats

**Coverage Matrix (`ofn_mat`):**

Tab-delimited file with header. For each library, two columns are generated:
- `cov_N`: Per-base coverage (number of reads / segment length)
- `var_N`: Variance of coverage

Values are formatted with 3 decimal places.

Example (4 libraries):
```
segment	cov_1	var_1	cov_2	var_2	cov_3	var_3	cov_4	var_4
seg_001	12.345	12.345	8.234	8.234	15.678	15.678	10.123	10.123
seg_002	6.789	6.789	5.432	5.432	7.890	7.890	6.543	6.543
seg_003	20.123	20.123	18.456	18.456	22.789	22.789	19.012	19.012
```

**Segment FASTA (`ofn_fasta`):**

FASTA file containing the extracted sequences for each segment based on the coordinates in the segment table.

Example:
```
>seg_001
ATCGATCGATCG...
>seg_002
GCTAGCTAGCTA...
>seg_003
TTAATTAATTAA...
```

## Algorithm

For each segment and library:
1. Load alignments from the library's ALN file
2. Query alignments intersecting the segment interval (1-based coordinates converted to 0-based half-open internally)
3. Count unique reads overlapping the segment
4. Calculate per-base coverage: `coverage = read_count / segment_length`

The command processes all segments for each library before moving to the next library to minimize file I/O.

## Use Cases

- **Prepare input for MetaBAT binning** with coverage information across multiple libraries
- **Generate segment-level coverage profiles** for metagenomic assemblies
- **Compare coverage patterns** across different conditions or timepoints
- **Extract segment sequences** along with their coverage statistics

## Notes

- Segment coordinates are 1-based and inclusive (standard biological convention)
- Coverage values represent per-base coverage (reads per base position)
- Progress messages are printed every 1000 segments processed
- The output maintains the same segment order as the input segment table

