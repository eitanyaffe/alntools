# alntools

`alntools` is a specialized toolkit for efficiently working with read alignments. It creates a compact binary representation of alignments (PAF format) and provides powerful querying capabilities to analyze specific genomic intervals. Key features include:

- **Fast binary storage** of read alignments from PAF format with mutation encoding
- **Three query modes** for flexible analysis:
  - **Full mode**: Retrieves complete read, alignment and mutation details with read-based height calculations for stacked visualization
  - **Pileup mode**: Provides position-by-position mutation summaries for variant analysis
  - **Bin mode**: Generates binned coverage statistics with segregating sites analysis and mutation rate categorization
- **Coverage analysis** for comprehensive alignment statistics and identification of unaligned regions
- **Break detection** for identifying positions with excessive read start/end clustering using statistical testing
- **R interface** for seamless integration with analysis workflows in R

This makes alntools useful for visualizing and investigating read coverage and mutation patterns across regions of interest.

## Installation

### Dependencies

*   A modern C++ compiler (supporting C++17)
*   gnumake
*   OpenMP library (for threading support in bin queries)
    *   **Linux**: Usually included with gcc (`libgomp`)
    *   **macOS**: Install via Homebrew: `brew install libomp`

#### GPU Acceleration (Optional)

*   **macOS with Apple Silicon**: Xcode Command Line Tools for Metal GPU acceleration
    *   Install: `xcode-select --install`
    *   Provides significant speedups for bin queries on large datasets
    *   Falls back to CPU automatically if not available

Tested on macOS 13.3.1 and Ubuntu 20.04.

#### For R Interface (Optional)

*   R (tested with 4.0+)
*   Rcpp package (`install.packages("Rcpp")`)

### Building

1.  Clone the repository:
    ```bash
    git clone https://github.com/eitanyaffe/alntools.git
    ```
2.  Compile the code:
    ```bash
    cd alntools
    make
    ```
    The executable will be located in `bin/macos/alntools` or `bin/linux/alntools`, depending on your system.
3.  (Optional) Install the executable:
    ```bash
    make install
    ```
    This copies the executable to /usr/local/bin.
4.  (Optional) Run tests:
    ```bash
    make test
    ```

## File Formats

For detailed information about all input and output file formats used by alntools, see the [File Formats documentation](FILE_FORMATS.md).

## CLI Commands

### 1. construct

Creates a binary `.aln` file from a PAF file, which stores alignment data efficiently for later queries.

```bash
alntools construct -ifn_paf <input.paf> -ofn <output.aln> [options]
```

**Mandatory Arguments:**
* `-ifn_paf <fn>`: Input alignment PAF file.
* `-ofn <fn>`: Path for the output ALN file.

**Optional Arguments:**
* `-verify <T|F>`: Verify PAF alignments against sequence files (default: `false`).
* `-ifn_reads <fn>`: Input read FASTQ file (required if `-verify T`).
* `-ifn_contigs <fn>`: Input contig FASTA file (required if `-verify T`).
* `-max_reads <int>`: Process only the first N alignments (0 means all, default: `0`).
* `-quit_on_error <T|F>`: Exit immediately if an error is encountered during parsing or verification (default: `true`).

**Example:**
```bash
# Basic construction without verification
mkdir output
alntools construct -ifn_paf examples/align_100.paf -ofn output/test.aln
```

### 2. info

Provides basic statistics about an ALN file.

```bash
alntools info -ifn <input.aln>
```

**Mandatory Arguments:**
* `-ifn <fn>`: Input ALN file.

**Example:**
```bash
alntools info -ifn output/test.aln
```

**Output Information:**
- Total alignments
- Total reads
- Average alignment length
- Total mutations
- Average mutations per alignment

### 3. query

Query the ALN file using different modes for specific contig intervals.

```bash
alntools query -ifn_aln <input.aln> -ifn_intervals <intervals.txt> -ofn_prefix <output_prefix> -mode <full|pileup|bin> [options]
```

**Mandatory Arguments:**
* `-ifn_aln <fn>`: Input ALN file.
* `-ifn_intervals <fn>`: Input tab-delimited file with query intervals (format: `contig start end`).
* `-ofn_prefix <fn>`: Output prefix for result files.
* `-mode <string>`: Query mode, one of:
  - `full`: Return detailed alignment and mutation data.
  - `pileup`: Return aggregated mutation data for positions.
  - `bin`: Return binned summaries of alignments.

**Optional Arguments (depending on mode):**
* `-pileup_mode <string>`: For pileup mode, options are:
  - `all`: Report all positions within query intervals.
  - `covered`: Report only positions with read coverage (default).
  - `mutated`: Report only positions with mutations.
* `-binsize <int>`: For bin mode, size of bins in bp (default: `100`).
* `-seg_threshold <double>`: For bin mode, threshold for segregating sites detection (default: `0.2`). Variants with frequency between this value and (1 - this value) are considered segregating.
* `-non_ref_threshold <double>`: For bin mode, threshold for non-reference sites detection (default: `0.9`). Variants with frequency above this value are considered non-reference.
* `-height_style <string>`: For full mode, how to calculate read height:
  - `by_coord_left`: Minimize overlap between reads, sort by start position (default).
  - `by_coord_right`: Minimize overlap between reads, sort by end position.
  - `by_mutations`: Arrange by mutation density.

**Alignment Filtering Options (all modes):**
* `-clip_mode <string>`: Clipping mode for filtering alignments:
  - `all`: Include all alignments (default).
  - `complete`: Only alignments covering the entire read.
  - `allow_one_side_clip`: Allow clipping on one side only.
  - `only_one_side_clipped`: Only alignments clipped on exactly one side.
  - `only_two_side_clipped`: Only alignments clipped on both sides.
  - `only_clipped`: Only alignments clipped on one or both sides.
* `-clip_margin <int>`: Margin in bases for clipping detection (default: `10`).
* `-min_mutations_percent <double>`: Minimum mutations percentage (default: `0.0`).
* `-max_mutations_percent <double>`: Maximum mutations percentage (default: `10.0`).
* `-min_alignment_length <int>`: Minimum alignment length in read coordinates (default: `0`).
* `-max_alignment_length <int>`: Maximum alignment length in read coordinates (default: `0`, no limit).
* `-min_indel_length <int>`: Minimum indel length to include in mutation density calculations (default: `3`). Shorter indels are filtered out to reduce noise from sequencing artifacts.

**Example of full query mode**

```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_large.txt \
   -ofn_prefix output/query -mode full
```

This produces three output files:
- `output/query_alignments.tsv`: Detailed alignment information with heights inherited from reads
- `output/query_mutations.tsv`: Mutation details for each alignment
- `output/query_reads.tsv`: Read statistics and height assignments (heights calculated on reads first)

**Example of bin query mode**

```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_small.txt \
   -ofn_prefix output/query -mode bin -binsize 1000 -seg_threshold 0.2 -non_ref_threshold 0.9 -num_threads 4
```

**Example with alignment filtering**

```bash
# Filter alignments by length and mutation rate, exclude short indels
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_small.txt \
   -ofn_prefix output/query_filtered -mode full \
   -min_alignment_length 1000 -max_alignment_length 50000 \
   -min_mutations_percent 0.1 -max_mutations_percent 5.0 \
   -min_indel_length 5 \
   -clip_mode complete
```

This produces a single output file:
- `output/query_bins.tsv`: Binned alignment statistics including:
  - Basic bin information: `contig`, `bin_start`, `bin_end`, `bin_length`
  - Coverage statistics: `sequenced_bp`, `mutation_count`
  - **Segregating sites analysis**: `seg_sites_density` (sites per bp)
  - **Non-reference sites analysis**: `non_ref_sites_density` (sites per bp)
  - **Mutation distance categories**: `dist_none` (0 mutations), `dist_5` (1e-5 to 1e-4 per bp), `dist_4` (1e-4 to 1e-3 per bp), `dist_3` (1e-3 to 1e-2 per bp), `dist_2` (1e-2 to 1e-1 per bp), `dist_1_plus` (>1e-1 per bp)

#### GPU Acceleration for Bin Queries

On macOS with Apple Silicon, bin queries can leverage GPU acceleration using Metal for significant performance improvements:

```bash
# GPU acceleration is automatically enabled when Metal is available
# The tool will automatically detect and use GPU when possible
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_large.txt \
   -ofn_prefix output/query -mode bin -binsize 1000
```

**Performance Benefits:**
- 1000x faster processing for large datasets (>10K alignments)
- Automatic fallback to CPU if GPU not available
- Bit-for-bit identical results to CPU implementation

**How GPU Acceleration Works:**
The bin query processes thousands of read alignments to calculate coverage and mutation statistics across genomic intervals. Traditionally, this requires iterating through each alignment sequentially on the CPU. GPUs excel at this type of work because they can process many alignments simultaneously using hundreds of parallel cores.

The implementation converts alignment data into a simplified format suitable for GPU processing, then uses Apple's Metal framework to dispatch the work across the GPU cores. Each GPU thread processes a subset of alignments, updating shared counters for sequenced base pairs, mutation counts, and distance categories using atomic operations to prevent conflicts. The results are then transferred back and combined with additional CPU-computed metrics (like segregating sites density) to ensure complete compatibility with the original CPU implementation.

**Example of pile-up query mode**

```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_small.txt \
   -ofn_prefix output/query -mode pileup -pileup_mode mutated
```

### 4. coverage

Generates comprehensive alignment coverage statistics for the entire ALN file, including detailed analysis of aligned and non-aligned regions.

```bash
alntools coverage -ifn <input.aln> -ofn_prefix <output_prefix>
```

**Mandatory Arguments:**
* `-ifn <fn>`: Input ALN file.
* `-ofn_prefix <fn>`: Output prefix for result files.

**Example:**
```bash
alntools coverage -ifn output/test.aln -ofn_prefix output/coverage_stats
```

**Output Files:**
- `{prefix}_coverage.txt`: Main statistics table with tab-delimited format containing:
  - `contig_count`: Number of contigs
  - `read_count`: Number of reads  
  - `alignment_count`: Number of alignments
  - `total_read_bp`: Total read base pairs
  - `total_assembly_bp`: Total assembly base pairs
  - `contig_n50`: N50 statistic of contig lengths (length of shortest contig such that 50% of assembly is in contigs of that length or longer)
  - `aligned_count`: Number of reads with at least one alignment
  - `non_aligned_read_bp`: Total base pairs in unaligned read regions
  - `non_aligned_contig_bp`: Total base pairs in unaligned contig regions
- `{prefix}_non_aligned_reads.txt`: Non-aligned read intervals with columns: read_id, start, end, length (1-based coordinates)
- `{prefix}_non_aligned_contigs.txt`: Non-aligned contig intervals with columns: contig_id, start, end, length (1-based coordinates)

**Use Cases:**
- Assess overall alignment quality and coverage
- Identify regions with poor alignment coverage
- Calculate genome-wide alignment statistics
- Find specific genomic intervals that lack read coverage

### 5. breaks

Identifies positions with an excessive number of reads that start or end exactly at that position, which may indicate structural variations, assembly artifacts, or other genomic features causing read clustering.

```bash
alntools breaks -ifn <input.aln> -ofn <output.tsv> -window <int> -pval <double> [options]
```

**Mandatory Arguments:**
* `-ifn <fn>`: Input ALN file.
* `-ofn <fn>`: Output TSV file for break positions.
* `-window <int>`: Window size for background calculation (e.g., 1000).
* `-pval <double>`: P-value threshold for reporting significant positions (e.g., 0.05).

**Optional Arguments:**
* `-min_reads <int>`: Minimum number of reads required to test a position (default: 1).

**Examples:**
```bash
# Basic usage with default min_reads (1)
alntools breaks -ifn output/test.aln -ofn output/breaks.tsv -window 1000 -pval 0.05

# Only test positions with at least 5 supporting reads
alntools breaks -ifn output/test.aln -ofn output/breaks.tsv -window 1000 -pval 0.05 -min_reads 5
```

**Output Format:**
The output TSV file contains the following columns (sorted by contig, then position):
- `contig`: Contig identifier
- `position`: 1-based position on the contig
- `orientation`: Either "left" (read starts) or "right" (read ends)
- `t`: Number of reads starting/ending at this position
- `e`: Expected number of reads per position (background rate)
- `enrichment`: Observed/expected ratio (t/e)
- `pval`: Raw p-value from binomial test
- `qval`: Benjamini-Hochberg corrected q-value

**Algorithm:**
- Only counts alignments where the contig boundary corresponds to the actual read start (`read_start == 0`) or read end (`read_end == read_length`)
- Filters positions to only test those with at least `min_reads` supporting reads
- For each position, compares observed counts to expected counts in a sliding window using a binomial test
- Calculates enrichment as observed/expected ratio (t/e)
- Applies Benjamini-Hochberg correction for multiple testing
- Reports positions with q-value ≤ threshold, sorted by contig and position

## R Interface

`alntools` provides an R interface for constructing, loading, and querying alignment stores.

### Loading the R Interface

**CPU-only version:**
```R
library(Rcpp)
sourceCpp("cpp/aln_cpu_R.cpp")
```

**GPU-enabled version (Apple Silicon):**
```R
library(Rcpp)
# Set environment for GPU compilation
Sys.setenv(PKG_CPPFLAGS = "-I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp -DMETAL_SUPPORT")
Sys.setenv(PKG_LIBS = "-L/opt/homebrew/opt/libomp/lib -lomp -framework Metal -framework MetalPerformanceShaders -framework Foundation")
Sys.setenv(PKG_CXXFLAGS = "-x objective-c++")
sourceCpp("cpp/aln_gpu_R.cpp")
```

### Available Functions

#### 1. Construction and Loading

```R
# input PAF file
paf_file <- "examples/align_100_dense.paf"

# output ALN file
aln_file <- "output/dense.aln"

# Construct alignment store from PAF file
aln <- aln_construct(paf_file, max_reads = 0)

# Save alignment store to file
aln_save(aln, aln_file)

# Load existing alignment store
aln <- aln_load(aln_file)
```

#### 2. Querying

```R
# input interval file
interval_file <- "examples/intervals_dense.txt"
intervals <- read.delim(interval_file)

# Bin query with threading support and optional GPU acceleration
bin_results <- aln_query_bin(aln, intervals, binsize, seg_threshold = 0.2, non_ref_threshold = 0.9, num_threads = 0, 
                            clip_mode_str = "all", clip_margin = 10, min_mutations_percent = 0.0, max_mutations_percent = 10.0,
                            min_alignment_length = 0, max_alignment_length = 0, min_indel_length = 3, use_gpu = FALSE)

# report_mode options: "all", "covered", "mutated"
report_mode <- "all"

# Pileup query
pileup_results <- aln_query_pileup(aln, intervals, report_mode, clip_mode_str = "all", clip_margin = 10, 
                                  min_mutations_percent = 0.0, max_mutations_percent = 10.0,
                                  min_alignment_length = 0, max_alignment_length = 0, min_indel_length = 3)

# height_style options: "by_coord_left", "by_coord_right", "by_mutations"
height_style <- "by_coord_left"
# Full query
full_results <- aln_query_full(aln, intervals, height_style, max_alignments = 0, clip_mode_str = "all", clip_margin = 10,
                              min_mutations_percent = 0.0, max_mutations_percent = 10.0,
                              min_alignment_length = 0, max_alignment_length = 0, min_indel_length = 3)
# Returns a list with $alignments, $mutations, and $reads dataframes

# Break position detection
breaks_results <- aln_find_breaks(aln, window_size = 1000, p_threshold = 0.05, min_reads = 1)
# Returns a dataframe with columns: contig, position, orientation, t, e, enrichment, pval, qval
# Sorted by contig, then position
```

### Example R Script

```R
# Load required packages
library(Rcpp)

# Compile and load C++ code
sourceCpp("cpp/aln_R.cpp")

# Parameters
paf_file <- "examples/align_100_dense.paf"
intervals_file <- "examples/intervals_dense.txt"
binsize <- 10
output_prefix <- "output/results"

# Construct and save alignment
aln_file <- paste0(output_prefix, ".aln")
aln <- aln_construct(paf_file)
aln_save(aln, aln_file)

# Load from file
aln <- aln_load(aln_file)

# Load intervals
intervals <- read.table(intervals_file, header = TRUE)

# Run queries
bin_results <- aln_query_bin(aln, intervals, binsize, seg_threshold = 0.2, non_ref_threshold = 0.9, num_threads = 0,
                            min_alignment_length = 1000, min_indel_length = 5)  # Filter alignments < 1kb, short indels < 5bp
pileup_results <- aln_query_pileup(aln, intervals, "covered", max_alignment_length = 50000, min_indel_length = 3)  # Filter alignments > 50kb, short indels < 3bp
full_results <- aln_query_full(aln, intervals, "by_mutations", clip_mode_str = "complete", min_indel_length = 3)  # Only complete alignments, filter short indels
breaks_results <- aln_find_breaks(aln, window_size = 1000, p_threshold = 0.05, min_reads = 3)

# Save results
write.table(bin_results, file = paste0(output_prefix, "_bins.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pileup_results, file = paste0(output_prefix, "_pileup.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_results$alignments, file = paste0(output_prefix, "_alignments.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_results$mutations, file = paste0(output_prefix, "_mutations.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_results$reads, file = paste0(output_prefix, "_reads.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(breaks_results, file = paste0(output_prefix, "_breaks.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
```

## Running Tests

The repository includes tests that demonstrate functionality and verify correctness:

```bash
# Run all CLI tests 
make test

# Run specific test groups
make test_basic     # Basic functionality
make test_query_all # All query modes
make test_coverage  # Coverage analysis
make test_R_all     # R interface tests
```

For more details on test scenarios, see `test.mk`.
