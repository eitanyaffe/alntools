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
3.  (Optional) install the executable:
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
