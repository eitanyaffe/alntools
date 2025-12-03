# query

Query the ALN file using different modes for specific contig intervals.

## Syntax

```bash
alntools query -ifn_aln <input.aln> -ifn_intervals <intervals.txt> -odir <output_dir> -mode <full|pileup|bin|consensus|variants> [options]
```

## Parameters

**Mandatory Arguments:**
* `-ifn_aln <fn>`: Input ALN file.
* `-ifn_intervals <fn>`: Input tab-delimited file with query intervals (format: `contig start end`).
* `-odir <dir>`: Output directory for result files.
* `-mode <string>`: Query mode, one of:
  - `full`: Return detailed alignment and mutation data.
  - `pileup`: Return aggregated mutation data for positions.
  - `bin`: Return binned summaries of alignments.
  - `consensus`: Return high-frequency variants above a consensus threshold.
  - `variants`: Multi-library variant calling across multiple ALN files.

**Mode-Specific Options:**

**Pileup Mode:**
* `-pileup_mode <string>`: Options are:
  - `all`: Report all positions within query intervals.
  - `covered`: Report only positions with read coverage (default).
  - `mutated`: Report only positions with mutations.

**Bin Mode:**
* `-binsize <int>`: Size of bins in bp (default: `100`).
* `-seg_threshold <double>`: Threshold for segregating sites detection (default: `0.2`). Variants with frequency between this value and (1 - this value) are considered segregating.
* `-non_ref_threshold <double>`: Threshold for non-reference sites detection (default: `0.9`). Variants with frequency above this value are considered non-reference.

**Consensus Mode:**
* `-consensus_threshold <double>`: Frequency threshold for reporting variants (default: `0.9`). Only variants with frequency ≥ this value are reported.
* `-min_consensus_coverage <int>`: Minimum coverage required for reporting variants (default: `5`). Only variants with coverage ≥ this value are considered.

**Variants Mode:**
* `-ifn_libraries <fn>`: Input tab-delimited file with library definitions (format: `id fn`). **Required for variants mode**.
* `-min_variants_variant_support <int>`: Minimum variant support across all libraries (default: `3`). Total reads supporting the variant across all libraries.
* `-min_variants_library_support <int>`: Minimum number of libraries with the variant (default: `1`). Number of libraries that must contain the variant.
* `-min_variants_coverage_support <int>`: Minimum total coverage across all libraries (default: `10`). Total coverage at the variant position across all libraries.

**Full Mode:**
* `-height_style <string>`: How to calculate read height:
  - `by_coord_left`: Minimize overlap between reads, sort by start position (default).
  - `by_coord_right`: Minimize overlap between reads, sort by end position.
  - `by_mutations`: Arrange by mutation density.
* `-chunk_type <string>`: How to define chunks (stretches of alignments) for height calculation:
  - `read`: Entire read forms one chunk (equivalent to read-based heights).
  - `alignment`: Each alignment forms its own chunk (maximum granularity).
  - `break_on_overlap`: Start new chunk when alignments overlap (default).
  - `break_on_gap`: Start new chunk when gap between alignments exceeds max_margin.
* `-max_margin <int>`: Maximum margin tolerance for chunk detection in read coordinates (default: `10`).

**Alignment Filtering Options (all modes):**
* `-clip_mode <string>`: Clipping mode for filtering alignments:
  - `all`: Include all alignments (default).
  - `complete`: Only alignments covering the entire read.
  - `allow_one_side_clip`: Allow clipping on one side only.
  - `only_one_side_clipped`: Only alignments clipped on exactly one side.
  - `only_two_side_clipped`: Only alignments clipped on both sides.
  - `only_clipped`: Only alignments clipped on one or both sides.
  - `local_align`: Only locally aligned reads (alignments that are part of reads that aligned on both sides to the current contig).
* `-clip_margin <int>`: Margin in bases for clipping detection (default: `10`).
* `-min_mutations_percent <double>`: Minimum mutations percentage (default: `0.0`).
* `-max_mutations_percent <double>`: Maximum mutations percentage (default: `10.0`).
* `-min_alignment_length <int>`: Minimum alignment length in read coordinates (default: `0`).
* `-max_alignment_length <int>`: Maximum alignment length in read coordinates (default: `0`, no limit).
* `-min_indel_length <int>`: Minimum indel length to include in mutation density calculations (default: `3`). Shorter indels are filtered out to reduce noise from sequencing artifacts.

**Threading Options:**
* `-num_threads <int>`: For bin and consensus modes, number of threads to use (default: `0`, auto-detect).

## Examples

**Full Query Mode:**
```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_large.txt \
   -odir output/query -mode full
```

**Bin Query Mode:**
```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_small.txt \
   -odir output/query -mode bin -binsize 1000 -seg_threshold 0.2 -non_ref_threshold 0.9 -num_threads 4
```

**Pileup Query Mode:**
```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_small.txt \
   -odir output/query -mode pileup -pileup_mode mutated
```

**Consensus Query Mode:**
```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_small.txt \
   -odir output/query -mode consensus \
   -consensus_threshold 0.9 -min_consensus_coverage 5 -num_threads 4
```

**Variants Query Mode (Multi-library):**
```bash
# Create libraries table file
echo -e "id\tfn" > libraries.txt
echo -e "lib1\toutput/sample1.aln" >> libraries.txt
echo -e "lib2\toutput/sample2.aln" >> libraries.txt
echo -e "lib3\toutput/sample3.aln" >> libraries.txt

# Run variants query
alntools query -mode variants \
   -ifn_libraries libraries.txt \
   -ifn_intervals examples/intervals_small.txt \
   -odir output/query_variants \
   -min_variants_variant_support 5 \
   -min_variants_library_support 2 \
   -min_variants_coverage_support 20
```

**Alignment Filtering Example:**
```bash
# Filter alignments by length and mutation rate, exclude short indels
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_small.txt \
   -odir output/query_filtered -mode full \
   -min_alignment_length 1000 -max_alignment_length 50000 \
   -min_mutations_percent 0.1 -max_mutations_percent 5.0 \
   -min_indel_length 5 \
   -clip_mode complete
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command
- **Intervals file**: See [File Formats](../file-formats.md#intervals-format)
- **Libraries file**: Tab-delimited file with columns `id` and `fn` (for variants mode)

## Output File Formats

See [File Formats](../file-formats.md) for detailed specifications of all output formats:

**Full Mode Output (in odir):**
- `full_alignments.txt`: Detailed alignment information with heights
- `full_mutations.txt`: Mutation details for each alignment
- `full_chunks.txt`: Chunk statistics and height assignments

**Pileup Mode Output (in odir):**
- `pileup_table.txt`: Position-by-position mutation summaries

**Bin Mode Output (in odir):**
- `bin_table.txt`: Binned alignment statistics with segregating sites analysis

**Consensus Mode Output (in odir):**
- `consensus_table.txt`: High-frequency variant calls

**Variants Mode Output (in odir):**
- `variant_table.txt`: Main variant table with library counts
- `variant_support.txt`: Read support matrix (variants × libraries)
- `variant_coverage.txt`: Coverage matrix (variants × libraries)

## Gene Annotation

All query modes support optional gene annotation to classify variants as genic or intergenic:

**Gene Annotation Options:**
* `-ifn_gene_table <fn>`: Tab-delimited gene table with columns: `gene`, `contig`, `start`, `end`, `strand`, `desc`
* `-ifn_codon_table <fn>`: Codon table file or name (e.g., `table11` for standard genetic code)
* `-ifn_reference_fasta <fn>`: Reference FASTA file for proper codon analysis
* `-use_genes <T|F>`: Enable gene annotation (requires gene table, codon table, and reference FASTA)

**Example with Gene Annotation:**
```bash
alntools query -mode variants \
   -ifn_libraries libraries.txt \
   -ifn_intervals intervals.txt \
   -odir output/annotated_variants \
   -ifn_gene_table genes.txt \
   -ifn_codon_table table11 \
   -ifn_reference_fasta reference.fa \
   -use_genes T
```

**Additional Output Files (when gene annotation is enabled):**
- Main output files include an additional `is_genic` column (true/false)
- `variant_genic.txt`: Detailed annotation for variants within genes
- `variant_intergenic.txt`: Information for variants between genes

## Notes

- The query command is the primary analysis tool for ALN files
- Different modes are optimized for different types of analysis
- Use alignment filtering options to focus on high-quality alignments
- Threading is available for computationally intensive bin and consensus modes
- Gene annotation provides functional context for identified variants
