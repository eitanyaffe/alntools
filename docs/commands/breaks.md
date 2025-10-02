# breaks

Identifies positions with an excessive number of reads that start or end exactly at that position, which may indicate structural variations, assembly artifacts, or other genomic features causing read clustering.

## Syntax

```bash
alntools breaks -ifn <input.aln> -ofn <output.txt> -window <int> -pval <double> [options]
```

## Parameters

**Mandatory Arguments:**
* `-ifn <fn>`: Input ALN file.
* `-ofn <fn>`: Output TSV file for break positions.
* `-window <int>`: Window size for background calculation (e.g., 1000).
* `-pval <double>`: P-value threshold for reporting significant positions (e.g., 0.05).

**Optional Arguments:**
* `-min_reads <int>`: Minimum number of reads required to test a position (default: 1).

## Examples

```bash
# Basic usage with default min_reads (1)
alntools breaks -ifn output/test.aln -ofn output/breaks.txt -window 1000 -pval 0.05

# Only test positions with at least 5 supporting reads
alntools breaks -ifn output/test.aln -ofn output/breaks.txt -window 1000 -pval 0.05 -min_reads 5
```

## Input File Formats

- **ALN file**: Binary alignment file created by the `construct` command

## Output File Formats

**Break Positions Table (`output.txt`):**

Tab-delimited format with columns (sorted by contig, then position):
- `contig`: Contig identifier
- `position`: 1-based position on the contig
- `orientation`: Either "left" (read starts) or "right" (read ends)
- `t`: Number of reads starting/ending at this position
- `e`: Expected number of reads per position (background rate)
- `enrichment`: Observed/expected ratio (t/e)
- `pval`: Raw p-value from binomial test
- `qval`: Benjamini-Hochberg corrected q-value

## Example Output

```
contig	position	orientation	t	e	enrichment	pval	qval
chr1	1500000	left	25	8.5	2.94	0.001	0.015
chr1	1500001	right	22	8.5	2.59	0.003	0.025
chr2	2800000	left	18	7.2	2.50	0.008	0.040
```

## Algorithm

The break detection algorithm works as follows:

1. **Filter alignments**: Only counts alignments where the contig boundary corresponds to the actual read start (`read_start == 0`) or read end (`read_end == read_length`)

2. **Position filtering**: Only tests positions with at least `min_reads` supporting reads

3. **Background calculation**: For each position, calculates expected counts using a sliding window of size `window`

4. **Statistical testing**: Compares observed counts to expected counts using a binomial test

5. **Enrichment calculation**: Computes observed/expected ratio (t/e)

6. **Multiple testing correction**: Applies Benjamini-Hochberg correction for multiple testing

7. **Output filtering**: Reports positions with q-value ≤ threshold, sorted by contig and position

## Use Cases

- **Structural variation detection**: Identify potential insertion, deletion, or rearrangement breakpoints
- **Assembly artifact identification**: Find positions where assembly errors cause abnormal read clustering
- **Quality control**: Assess alignment quality by identifying regions with unusual read start/end patterns
- **Genomic feature discovery**: Locate positions that may correspond to repetitive elements or other features causing read clustering

## Notes

- The algorithm specifically looks for true read boundaries, not internal alignment breaks
- Use appropriate window sizes based on your data characteristics (larger windows for sparser data)
- The `min_reads` parameter helps reduce noise from low-coverage positions
- Statistical correction accounts for multiple testing across many genomic positions
- Results should be interpreted in the context of your specific biological system and sequencing technology
