# alntools

`alntools` represents read alignments to metagenomic assemblies and allows quick querying of specific intervals within contigs to fetch all reads or read summaries such as mutation pileups.

## Installation

### Dependencies

*   A modern C++ compiler (supporting C++17)
*   `make`
*   `zlib` library (e.g., `zlib1g-dev` on Debian/Ubuntu, `zlib-devel` on Fedora/CentOS, or install via Homebrew on macOS: `brew install zlib`)

Tested on:
*   macOS 13.3.1
*   Ubuntu 20.04

### Building

1.  Clone the repository:
    ```bash
    git clone <repository-url>
    cd alntools
    ```
2.  Compile the code:
    ```bash
    make
    ```
    The executable will be located in `bin/macos/alntools` or `bin/linux/alntools` depending on your system.
3.  (Optional) Run tests:
    ```bash
    make test
    ```

## Input Format (PAF)

`alntools` processes alignment data in the Pairwise mApping Format (PAF). The following columns are required:

1.  **Query sequence name** (string)
2.  **Query sequence length** (int)
3.  **Query start** (int, 0-based)
4.  **Query end** (int, 0-based)
5.  **Relative strand** (char, '+' or '-')
6.  **Target sequence name** (string) - *This corresponds to the contig ID.*
7.  **Target sequence length** (int)
8.  **Target start on original strand** (int, 0-based)
9.  **Target end on original strand** (int, 0-based)
10. **Number of residue matches** (int) - *Required by PAF spec, checked by alntools, but value not directly used.*
11. **Alignment block length** (int) - *Required by PAF spec, checked by alntools, but value not directly used.*
12. **Mapping quality** (int) - *Required by PAF spec, checked by alntools, but value not directly used.*

Additionally, the **`cs:Z:` tag** (difference string) **must be present** in one of the optional fields (column 13 onwards). This tag encodes the base differences between the query (read) and the target (contig) and is essential for mutation analysis.

### Example PAF Line

```text
read_298    150 10  145 +   contig_12   5000    1500    1635    120 135 60  tp:A:P  cm:i:10 NM:i:5  cs:Z::10*at:5+gg:20-c:15*cg:70
```
*(Note: This is a generic example; ensure your `cs` tag format matches what your aligner produces and what `alntools` expects)*

## Data Representation

*(No specific details requested)*

## Usage

### Constructing the Alignment Store (`.aln` file)

The `construct` command reads a PAF file and creates a binary `.aln` file, which stores the alignment data efficiently for later queries.

```bash
./bin/<os>/alntools construct -ifn_paf <input.paf> -ofn <output.aln> [options]
```

**Mandatory Arguments:**

*   `-ifn_paf <fn>`: Input alignment PAF file.
*   `-ofn <fn>`: Path for the output ALN file.

**Optional Arguments:**

*   `-verify <T|F>`: Verify PAF alignments against sequence files (default: `false`).
*   `-ifn_reads <fn>`: Input read FASTQ file (required if `-verify T`).
*   `-ifn_contigs <fn>`: Input contig FASTA file (required if `-verify T`).
*   `-max_reads <int>`: Process only the first N alignments (0 means all, default: `0`).
*   `-quit_on_error <T|F>`: Exit immediately if an error is encountered during parsing or verification (default: `true`).

---

*(Query commands will be added later)* 