# coverage

Generates comprehensive alignment coverage statistics for the entire ALN file, including detailed analysis of aligned and non-aligned regions.

## Syntax

```bash
alntools coverage -ifn <input.aln> -ofn_prefix <output_prefix>
```

## Parameters

**Mandatory Arguments:**
* `-ifn <fn>`: Input ALN file.
* `-ofn_prefix <fn>`: Output prefix for result files.

## Example

```bash
alntools coverage -ifn output/test.aln -ofn_prefix output/coverage_stats
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command

## Output File Formats

**Main Statistics Table (`{prefix}_coverage.txt`):**

Tab-delimited format containing:
- `contig_count`: Number of contigs
- `read_count`: Number of reads  
- `alignment_count`: Number of alignments
- `total_read_bp`: Total read base pairs
- `total_assembly_bp`: Total assembly base pairs
- `contig_n50`: N50 statistic of contig lengths (length of shortest contig such that 50% of assembly is in contigs of that length or longer)
- `aligned_count`: Number of reads with at least one alignment
- `non_aligned_read_bp`: Total base pairs in unaligned read regions
- `non_aligned_contig_bp`: Total base pairs in unaligned contig regions

**Non-aligned Regions Tables:**
- `{prefix}_non_aligned_reads.txt`: Non-aligned read intervals with columns: `read_id`, `start`, `end`, `length` (1-based coordinates)
- `{prefix}_non_aligned_contigs.txt`: Non-aligned contig intervals with columns: `contig_id`, `start`, `end`, `length` (1-based coordinates)

## Example Output

**Coverage Statistics:**
```
contig_count	1250
read_count	10000
alignment_count	12500
total_read_bp	50000000
total_assembly_bp	45000000
contig_n50	15000
aligned_count	9850
non_aligned_read_bp	2500000
non_aligned_contig_bp	1800000
```

**Non-aligned Reads:**
```
read_id	start	end	length
read_001	1	150	150
read_002	500	750	250
read_003	1200	1500	300
```

## Use Cases

- **Assess overall alignment quality** and coverage across the entire dataset
- **Identify regions with poor alignment coverage** that may need additional sequencing or different alignment parameters
- **Calculate genome-wide alignment statistics** for quality control
- **Find specific genomic intervals** that lack read coverage for targeted analysis
- **Evaluate completeness** of alignment datasets before downstream analysis

## Notes

- This command analyzes the entire ALN file, not specific intervals
- The N50 statistic helps assess the contiguity of the assembly
- Non-aligned regions may indicate areas of poor sequence quality, repetitive regions, or gaps in coverage
- Use the non-aligned region files to identify specific areas that may need attention
- This is particularly useful for quality control before running detailed queries
