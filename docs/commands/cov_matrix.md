# cov_matrix

Computes per-base coverage and variance statistics for genomic segments across multiple libraries, producing a coverage matrix suitable for downstream binning analysis (e.g., MetaBAT).

## Syntax

```bash
alntools cov_matrix -ifn_libraries <library_table> \
                    -ifn_segments <segment_table> \
                    -ifn_fasta <contigs.fasta> \
                    -ofn_mat <coverage_matrix.txt> \
                    -ofn_fasta <segments.fasta> \
                    [-actual_nts <T|F>] \
                    [-should_create_fasta <T|F>] \
                    [-min_segment_length <int>] \
                    [-clip_mode <mode>] \
                    [-clip_margin <int>] \
                    [-min_mutations_percent <double>] \
                    [-max_mutations_percent <double>] \
                    [-min_alignment_length <int>] \
                    [-max_alignment_length <int>] \
                    [-min_indel_length <int>] \
                    [-variance_mode <mode>]
```

## Parameters

**Mandatory Arguments:**
* `-ifn_libraries <fn>`: Library table with `lib_id` and `aln_fn` columns.
* `-ifn_segments <fn>`: Segment table with `segment`, `contig`, `start`, `end` columns.
* `-ifn_fasta <fn>`: Contig fasta file.
* `-ofn_mat <fn>`: Output coverage matrix file.
* `-ofn_fasta <fn>`: Output segment fasta file.

**Optional Arguments - Output Control:**
* `-actual_nts <boolean>`: When `T` (default), write actual segment sequences to the fasta. When `F`, random nucleotides are emitted instead.
* `-should_create_fasta <boolean>`: When `T` (default), generate the segment fasta file. Set to `F` to skip fasta generation entirely.

**Optional Arguments - Segment Filtering:**
* `-min_segment_length <int>`: Minimum segment length for filtering (default: 1000). Segments shorter than this threshold will be excluded from both output files.

**Optional Arguments - Alignment Filtering:**
* `-clip_mode <mode>`: Clipping mode for alignment filtering (default: end_unique). Options:
  * `all`: Allow all alignments
  * `complete`: Alignment must cover all read from start to end
  * `allow_one_side_clip`: Allow clipped on one side (start at read start or read end)
  * `only_one_side_clipped`: Show only alignments clipped on one side
  * `only_two_side_clipped`: Show only alignments clipped on both sides
  * `only_clipped`: Show alignments clipped on one or both sides
  * `local_align`: Show only locally aligned reads (first/last alignments on same contig)
  * `end_unique`: Like `complete`, but also accepts an alignment clipped on one read side only, where the clip lies at the contig end (read overhanging the contig), provided no other alignment of the same read overlaps more than half of it in read coordinates (excludes secondary placements). Clips inside the contig and alignments clipped on both read sides are rejected
* `-clip_margin <int>`: Clipping margin in bases (default: 10). Used to determine if alignment covers read start/end.
* `-min_mutations_percent <double>`: Minimum mutations percentage (default: 0.0). Alignments with fewer mutations are excluded.
* `-max_mutations_percent <double>`: Maximum mutations percentage (default: 0.1). Alignments with more mutations are excluded.
* `-min_alignment_length <int>`: Minimum alignment length in read coordinates (default: 1000). Shorter alignments are excluded.
* `-max_alignment_length <int>`: Maximum alignment length in read coordinates (default: 0, no limit). Longer alignments are excluded when > 0.
* `-min_indel_length <int>`: Minimum indel length to include in mutation density calculations (default: 3). Smaller indels are excluded from mutation counts.

**Optional Arguments - Variance:**
* `-variance_mode <mode>`: How `var_N` is computed (default: empirical). Options:
  * `empirical`: Population variance of the per-base depth along the segment, computed from the filtered alignments by sweeping their start and end positions (zero-depth bases included). This is the per-base mean and variance MetaBAT2 expects as input
  * `poisson`: Variance equals the mean depth

## Examples

**Basic usage with default filtering:**

```bash
alntools cov_matrix -ifn_libraries libraries.txt \
                    -ifn_segments segments.txt \
                    -ifn_fasta contigs.fa \
                    -ofn_mat coverage_matrix.txt \
                    -ofn_fasta segments.fa
```

**Custom filtering for high-quality alignments:**

```bash
alntools cov_matrix -ifn_libraries libraries.txt \
                    -ifn_segments segments.txt \
                    -ifn_fasta contigs.fa \
                    -ofn_mat coverage_matrix.txt \
                    -ofn_fasta segments.fa \
                    -max_mutations_percent 0.05 \
                    -min_alignment_length 2000
```

**Strict end-to-end alignments with Poisson variance (behavior of versions before 1.02):**

```bash
alntools cov_matrix -ifn_libraries libraries.txt \
                    -ifn_segments segments.txt \
                    -ifn_fasta contigs.fa \
                    -ofn_mat coverage_matrix.txt \
                    -ofn_fasta segments.fa \
                    -clip_mode complete \
                    -variance_mode poisson
```

**Include all alignments without filtering:**

```bash
alntools cov_matrix -ifn_libraries libraries.txt \
                    -ifn_segments segments.txt \
                    -ifn_fasta contigs.fa \
                    -ofn_mat coverage_matrix.txt \
                    -ofn_fasta segments.fa \
                    -clip_mode all \
                    -min_alignment_length 0
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
- `cov_N`: Mean per-base depth (aligned bp within the segment / segment length)
- `var_N`: Per-base depth variance; with `-variance_mode empirical` measured along the segment, with `-variance_mode poisson` equal to `cov_N`

Values are formatted with 3 decimal places.

Example (4 libraries, empirical variance):
```
segment	cov_1	var_1	cov_2	var_2	cov_3	var_3	cov_4	var_4
seg_001	12.345	0.812	8.234	0.451	15.678	1.207	10.123	0.633
seg_002	6.789	7.214	5.432	6.018	7.890	8.335	6.543	7.102
seg_003	20.123	22.457	18.456	19.873	22.789	25.110	19.012	21.364
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
2. Initialize short indel counting for mutation density calculations
3. Initialize the read-to-alignments index if the clip mode needs it (`local_align`, `end_unique`)
4. Query alignments intersecting the segment interval (1-based coordinates converted to 0-based half-open internally)
5. Filter alignments based on specified criteria (clip mode, mutation percent, alignment length)
6. Sum the bases of each passing alignment that fall inside the segment
7. Mean depth: `cov = aligned_bases / segment_length`
8. Variance: with `empirical`, sort the start (+1) and end (-1) positions of the passing alignments within the segment and sweep them, accumulating `width * depth^2` over each stretch of constant depth; `var = sum / segment_length - cov^2`. With `poisson`, `var = cov`

The command processes all segments for each library before moving to the next library to minimize file I/O.

## Use Cases

- **Prepare input for MetaBAT binning** with coverage information across multiple libraries
- **Generate segment-level coverage profiles** for metagenomic assemblies
- **Compare coverage patterns** across different conditions or timepoints
- **Extract segment sequences** along with their coverage statistics

## Notes

- Segment coordinates are 1-based and inclusive (standard biological convention)
- Coverage values are mean per-base depth (aligned bases per segment base)
- With long reads, the empirical variance of a segment much shorter than the reads is small (reads span it), and MetaBAT2 floors variances below 1 at 1
- Progress messages are printed every 1000 segments processed
- The output maintains the same segment order as the input segment table
- Segments shorter than `min_segment_length` are filtered out during loading and excluded from both output files
- Filtering statistics are reported showing how many segments were excluded
- Alignment filtering uses the shared alignment filter (also used by `csegment_coverage`, `get_read_ids` and the query modes; `end_unique` is not accepted by `query`)
- Default alignment filters (clip_mode=end_unique, max_mutations_percent=0.1, min_alignment_length=1000) count end-to-end alignments plus uniquely placed reads overhanging a contig end, which avoids the depth deficit that `complete` causes within about one read length of contig ends
- Use `-clip_mode all -max_mutations_percent -1 -min_alignment_length 0` to disable all alignment filtering and include all alignments

