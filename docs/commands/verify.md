# verify

Verifies ALN file consistency by checking that alignments match the original read and contig sequences.

## Syntax

```bash
alntools verify -ifn_aln <input.aln> -ifn_reads <reads.fastq> -ifn_contigs <contigs.fasta> [options]
```

## Parameters

**Mandatory Arguments:**
* `-ifn_aln <fn>`: Input ALN file.
* `-ifn_reads <fn>`: Input reads FASTQ file.
* `-ifn_contigs <fn>`: Input contigs FASTA file.

**Optional Arguments:**
* `-max_reads <int>`: Verify only the first N alignments (0 means all, default: `0`).
* `-ofn_reads <fn>`: Output FASTQ file with reads used in alignments.
* `-ofn_contigs <fn>`: Output FASTA file with contigs used in alignments.

## Example

```bash
# Basic verification
alntools verify -ifn_aln output/test.aln \
  -ifn_reads examples/reads_100.fq \
  -ifn_contigs examples/contigs_100.fa

# Verify only first 100 alignments and save filtered sequences
alntools verify -ifn_aln output/test.aln \
  -ifn_reads examples/reads_100.fq \
  -ifn_contigs examples/contigs_100.fa \
  -max_reads 100 \
  -ofn_reads output/filtered_reads.fq \
  -ofn_contigs output/filtered_contigs.fa
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command
- **FASTQ file**: Standard FASTQ format containing the original read sequences
- **FASTA file**: Standard FASTA format containing the original contig sequences

## Output File Formats

- **Console output**: Verification results and statistics
- **FASTQ file**: Filtered reads (if `-ofn_reads` specified)
- **FASTA file**: Filtered contigs (if `-ofn_contigs` specified)

## Verification Process

The verification process checks consistency between the ALN file and original sequences by:

1. **Extracting alignment information** from the ALN file
2. **Applying mutations** from the PAF cs:Z: tags to contig segments
3. **Comparing mutated contig segments** with corresponding read segments
4. **Reporting mismatches** with detailed position information

## Example Output

```
Reading alignment file: output/test.aln
Processing alignment 1/1250...
Alignment is good.
Processing alignment 2/1250...
Mismatch found, fragment coordinate=45
read        : ATCGATCGATCG
contig_mut  : ATCGATCAATCG
contig_orig : ATCGATCAATCG
...
Verification complete. Total alignments processed: 1250
Bad alignments found: 3 out of 1250 (0.24%)
```

## Notes

- This command is essential for validating the correctness of ALN files
- Verification ensures that the mutation encoding in the ALN file accurately represents the differences between reads and contigs
- Use this command when debugging alignment issues or validating new PAF data
- The optional output files contain only the sequences that were actually used in the alignments
- Verification can be time-consuming for large datasets; use `-max_reads` for testing
