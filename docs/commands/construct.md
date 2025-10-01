# construct

Creates a binary `.aln` file from a PAF file, which stores alignment data efficiently for later queries.

## Syntax

```bash
alntools construct -ifn_paf <input.paf> -ofn <output.aln> [options]
```

## Parameters

**Mandatory Arguments:**
* `-ifn_paf <fn>`: Input alignment PAF file.
* `-ofn <fn>`: Path for the output ALN file.

**Optional Arguments:**
* `-verify <T|F>`: Verify PAF alignments against sequence files (default: `false`).
* `-ifn_reads <fn>`: Input read FASTQ file (required if `-verify T`).
* `-ifn_contigs <fn>`: Input contig FASTA file (required if `-verify T`).
* `-max_reads <int>`: Process only the first N alignments (0 means all, default: `0`).
* `-quit_on_error <T|F>`: Exit immediately if an error is encountered during parsing or verification (default: `true`).

## Example

```bash
# Basic construction without verification
mkdir output
alntools construct -ifn_paf examples/align_100.paf -ofn output/test.aln

# Construction with verification
alntools construct -ifn_paf examples/align_100.paf -ofn output/test.aln \
  -verify T -ifn_reads examples/reads_100.fq -ifn_contigs examples/contigs_100.fa

# Process only first 1000 alignments
alntools construct -ifn_paf examples/align_100.paf -ofn output/test.aln -max_reads 1000
```

## Input File Formats

- **PAF file**: See [File Formats](../file-formats.md#paf-format) for detailed PAF format specification
- **FASTQ file**: Standard FASTQ format for reads (required only if verification is enabled)
- **FASTA file**: Standard FASTA format for contigs (required only if verification is enabled)

## Output File Formats

- **ALN file**: Binary format optimized for fast querying. Use the `info` command to view basic statistics or `extract` command to export to readable tables.

## Notes

- The PAF file must contain the `cs:Z:` tag (difference string) for mutation analysis
- Verification compares the PAF alignments against the original sequence files to ensure consistency
- The binary ALN format provides significant performance improvements for subsequent query operations
- Use `max_reads` parameter for testing with large datasets
