# alntools

`alntools` is a specialized toolkit for efficiently working with read alignments. It creates a compact binary representation of alignments (PAF format) and provides powerful querying capabilities to analyze specific genomic intervals. Key features include:

- **Fast binary storage** of read alignments from PAF format with mutation encoding
- **Genome rearrangement detection** for identifying large insertions, deletions, and inversions from read alignment patterns
- **Five query modes** for flexible analysis:
  - **Full mode**: Retrieves complete read, alignment and mutation details with read-based height calculations for stacked visualization
  - **Pileup mode**: Provides position-by-position mutation summaries for variant analysis
  - **Bin mode**: Generates binned coverage statistics with segregating sites analysis and mutation rate categorization
  - **Consensus mode**: Identifies high-frequency variants above a consensus threshold for variant calling
  - **Variants mode**: Multi-library variant calling across multiple ALN files with comprehensive variant filtering
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

## Commands

| Command | Description | Documentation |
|---------|-------------|---------------|
| [construct](docs/commands/construct.md) | Create binary ALN file from PAF file | [→](docs/commands/construct.md) |
| [info](docs/commands/info.md) | Show basic statistics for ALN file | [→](docs/commands/info.md) |
| [extract](docs/commands/extract.md) | Export ALN file to tab-delimited tables | [→](docs/commands/extract.md) |
| [verify](docs/commands/verify.md) | Verify ALN file consistency against original sequences | [→](docs/commands/verify.md) |
| [query](docs/commands/query.md) | Query ALN file with multiple analysis modes | [→](docs/commands/query.md) |
| [coverage](docs/commands/coverage.md) | Generate alignment coverage statistics | [→](docs/commands/coverage.md) |
| [breaks](docs/commands/breaks.md) | Find positions with excessive read start/end clustering | [→](docs/commands/breaks.md) |
| [rearrange](docs/commands/rearrange.md) | Detect genome rearrangements from read alignments | [→](docs/commands/rearrange.md) |

## Documentation

- **[File Formats](docs/file-formats.md)**: Detailed specifications for all input and output file formats
- **[R Interface](docs/r-interface.md)**: R interface documentation for integrating alntools into R workflows

## Quick Start

```bash
# 1. Create binary alignment file from PAF
alntools construct -ifn_paf examples/align_100.paf -ofn output/test.aln

# 2. View basic statistics
alntools info -ifn output/test.aln

# 3. Query specific intervals
alntools query -ifn_aln output/test.aln \
   -ifn_intervals examples/intervals_large.txt \
   -ofn_prefix output/query -mode full

# 4. Generate coverage statistics
alntools coverage -ifn output/test.aln -ofn_prefix output/coverage

# 5. Detect structural variations
alntools rearrange -ifn_aln output/test.aln -ofn_prefix output/rearrangements
```

## Running Tests

The repository includes tests that demonstrate functionality and verify correctness:

```bash
# Run all CLI tests 
make test

# Run specific test groups
make test_basic           # Basic functionality
make test_query_all       # All query modes (full, pileup, bin, consensus, variants)
make test_query_consensus # Consensus mode specifically
make test_query_variants  # Variants mode specifically
make test_coverage        # Coverage analysis
make test_genes           # Gene annotation functionality
make test_R_commands      # R interface tests
```

For more details on test scenarios, see `test.mk`.

## License

See [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.