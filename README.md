# alntools

`alntools` is a specialized toolkit for efficiently working with read alignments. It creates a compact binary representation of alignments (PAF format) and provides powerful querying capabilities to analyze specific genomic intervals. Key features include:

- **Fast binary storage** of read alignments from PAF format with mutation encoding
- **Three query modes** for flexible analysis:
  - **Full mode**: Retrieves complete alignment and mutation details with height calculations for stacked visualization
  - **Pileup mode**: Provides position-by-position mutation summaries for variant analysis
  - **Bin mode**: Generates binned coverage statistics for detecting coverage patterns
- **R interface** for seamless integration with analysis workflows in R

This makes alntools useful for visualizing and investigating read coverage and mutation patterns across regions of interest.

## Installation

### Dependencies

*   A modern C++ compiler (supporting C++17)
*   gnumake

Tested on macOS 13.3.1 and Ubuntu 20.04.

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
* `-height_style <string>`: For full mode, how to calculate alignment height:
  - `by_coord`: Minimize overlap between alignments (default).
  - `by_mutations`: Arrange by mutation density.

**Example of full query mode**

```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_large.txt \
   -ofn_prefix output/query -mode full
```

**Example of bin query mode**

```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_small.txt \
   -ofn_prefix output/query -mode bin -binsize 1000
```

**Example of pile-up query mode**

```bash
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_small.txt \
   -ofn_prefix output/query -mode pileup -pileup_mode mutated
```

## R Interface

`alntools` provides an R interface for constructing, loading, and querying alignment stores.

### Loading the R Interface

```R
library(Rcpp)
sourceCpp("cpp/aln_R.cpp")
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
# Bin query
bin_results <- aln_query_bin(aln, intervals, binsize)

# Pileup query
pileup_results <- aln_query_pileup(aln, intervals, report_mode)
# report_mode options: "all", "covered", "mutated"

# Full query
full_results <- aln_query_full(aln, intervals, height_style)
# height_style options: "by_coord", "by_mutations"
# Returns a list with $alignments and $mutations dataframes
```

#### 3. Visualization

```R
# Plot alignments
source("R/plot.r")
plot_alignments(alignments_df, output_file)
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
bin_results <- aln_query_bin(aln, intervals, binsize)
pileup_results <- aln_query_pileup(aln, intervals, "covered")
full_results <- aln_query_full(aln, intervals, "by_mutations")

# Save results
write.table(bin_results, file = paste0(output_prefix, "_bins.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pileup_results, file = paste0(output_prefix, "_pileup.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_results$alignments, file = paste0(output_prefix, "_alignments.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_results$mutations, file = paste0(output_prefix, "_mutations.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
```

## Running Tests

The repository includes tests that demonstrate functionality and verify correctness:

```bash
# Run all tests
make test

# Run specific test groups
make -f test.mk test_basic     # Basic functionality
make -f test.mk test_query_all # All query modes
make -f test.mk test_R_all     # R interface tests
```

For more details on test scenarios, see `test.mk`.
