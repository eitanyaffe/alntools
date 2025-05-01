# alntools

`alntools` represents read alignments to metagenomic assemblies and allows quick querying of specific intervals within contigs to fetch all reads or read summaries such as mutation pileups.

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

## Input Format (PAF)

`alntools` processes alignment data in the Pairwise mApping Format (PAF). The following columns are required:

| Column Index | Description                                      | Type  |
|--------------|--------------------------------------------------|-------|
| 1            | Query sequence name                              | string|
| 2            | Query sequence length                            | int   |
| 3            | Query start (0-based)                            | int   |
| 4            | Query end (0-based)                              | int   |
| 5            | Relative strand ('+' or '-')                     | char  |
| 6            | Target sequence name (contig ID)                 | string|
| 7            | Target sequence length                           | int   |
| 8            | Target start on original strand (0-based)        | int   |
| 9            | Target end on original strand (0-based)          | int   |

Additionally, the **`cs:Z:` tag** (difference string) **must be present** in one of the optional fields (column 13 onwards). This tag encodes the base differences between the query (read) and the target (contig) and is essential for mutation analysis.

### Example PAF Line

```text
read_298    150 10  145 +   contig_12   5000    1500    1635    120 135 60  cs:Z::10*at:5+gg:20-c:15*cg:70
```

For a complete example of a PAF file, see `examples/align_100.paf` in the repository.

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
alntools construct -ifn_paf examples/align_100.paf -ofn output/test.aln

# Construction with verification
alntools construct -ifn_paf examples/align_100.paf -ifn_reads examples/reads_100.fq -ifn_contigs examples/contigs_100.fa -verify T -ofn output/test_full.aln
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

**Examples:**
```bash
# Query in full mode
alntools query -ifn_aln output/test.aln -ifn_intervals examples/intervals_large.txt -ofn_prefix output/query -mode full

# Query in bin mode
alntools query -ifn_aln output/test.aln -ifn_intervals examples/intervals_small.txt -ofn_prefix output/query -mode bin -binsize 1000

# Query in pileup mode
alntools query -ifn_aln output/test.aln -ifn_intervals examples/intervals_small.txt -ofn_prefix output/query -mode pileup -pileup_mode mutated
```

### Intervals File Format

The intervals file is a tab-delimited file specifying regions to query:

| Column  | Description                           | Type  |
|---------|---------------------------------------|-------|
| contig  | Name of the contig                    | string|
| start   | Start position (1-based, inclusive)   | int   |
| end     | End position (1-based, inclusive)     | int   |

Example intervals file:
```
contig  start   end
ctg25860    120448  130902
ctg26175    230065  252386
```

## Query Output Formats

### 1. Full Mode Output

Produces two tab-delimited files:

**1. *_alignments.tsv:**

| Column         | Description                              | Type    |
|----------------|------------------------------------------|---------|
| alignment_index| Unique index for the alignment           | int     |
| read_id        | ID of the read                           | string  |
| contig_id      | ID of the contig                         | string  |
| read_start     | Start position on read (0-based)         | int     |
| read_end       | End position on read (0-based)           | int     |
| contig_start   | Start position on contig (0-based)       | int     |
| contig_end     | End position on contig (0-based)         | int     |
| is_reverse     | Whether alignment is on reverse strand   | boolean |
| cs_tag         | CIGAR string encoding differences        | string  |
| num_mutations  | Number of mutations in the alignment     | int     |
| height         | Vertical position for visualization      | int     |

**2. *_mutations.tsv:**

| Column         | Description                              | Type    |
|----------------|------------------------------------------|---------|
| alignment_index| Index of parent alignment                | int     |
| read_id        | ID of the read                           | string  |
| contig_id      | ID of the contig                         | string  |
| type           | Mutation type (SUB/INS/DEL)              | string  |
| position       | Position on contig (0-based)             | int     |
| desc           | Description of the mutation              | string  |
| height         | Vertical position for visualization      | int     |

### 2. Pileup Mode Output

Produces *_pileup.tsv:

| Column   | Description                                 | Type   |
|----------|---------------------------------------------|--------|
| contig   | Contig ID                                   | string |
| position | Position on contig (1-based)                | int    |
| variant  | Variant type (REF or mutation description)  | string |
| count    | Count of occurrences                        | int    |
| coverage | Total read coverage at this position        | int    |
| cumsum   | Cumulative count                            | int    |

### 3. Bin Mode Output

Produces *_bins.tsv:

| Column          | Description                               | Type   |
|-----------------|-------------------------------------------|--------|
| contig          | Contig ID                                 | string |
| start           | Bin start position                        | int    |
| end             | Bin end position                          | int    |
| length          | Bin length                                | int    |
| read_count      | Number of alignments in bin               | int    |
| mutation_count  | Number of mutations in bin                | int    |

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
# Construct alignment store from PAF file
aln <- aln_construct(paf_file, max_reads = 0)

# Save alignment store to file
aln_save(aln, output_file)

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
paf_file <- "examples/align_100.paf"
intervals_file <- "examples/intervals_small.txt"
binsize <- 1000
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
