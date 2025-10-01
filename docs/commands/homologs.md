# homologs

Performs kmer-based homology search to find regions in an assembly that share sequence similarity with a query region. This command extracts kmers from a specified query region and searches the entire assembly to identify homologous regions based on kmer coverage thresholds.

## Syntax

```bash
alntools homologs -ifn_fasta <assembly.fasta> -contig <query_contig> -start <start_pos> -end <end_pos> -ofn <output.txt> [options]
```

## Parameters

**Mandatory Arguments:**
* `-ifn_fasta <fn>`: Input assembly FASTA file.
* `-contig <string>`: Query contig name.
* `-start <int>`: Query start position (1-based).
* `-end <int>`: Query end position (1-based, inclusive).
* `-ofn <fn>`: Output regions file.

**Optional Arguments:**
* `-k <int>`: Kmer size (default: 21).
* `-num_kmers <int>`: Number of kmers to extract from query region (default: 10).
* `-threshold <double>`: Minimum percentage of kmers required for a match (default: 80.0).

## Example

```bash
# Find regions similar to positions 5000-15000 on contig_A
alntools homologs \
  -ifn_fasta genome.fasta \
  -contig contig_A \
  -start 5000 \
  -end 15000 \
  -k 25 \
  -num_kmers 20 \
  -threshold 75 \
  -ofn homologous_regions.txt
```

## Input File Formats

- **FASTA file**: Standard multi-sequence FASTA format containing the assembly sequences

## Output File Formats

**Regions Table (`output.txt`):**

Tab-delimited format compatible with mview regions, containing:
- `assembly`: Assembly identifier (always "assembly")
- `contig`: Contig name containing the homologous region
- `start`: Start position of the homologous region (1-based)
- `end`: End position of the homologous region (1-based, inclusive)
- `desc`: Description including region length and kmer coverage percentage
- `id`: Unique identifier for the homologous region (homolog_1, homolog_2, etc.)

## Example Output

**Homologous Regions:**
```
assembly	contig	start	end	desc	id
assembly	contig_B	12450	22890	kmer_match_length_10440_coverage_85pct	homolog_1
assembly	contig_C	5200	14100	kmer_match_length_8900_coverage_80pct	homolog_2
assembly	contig_D	18000	25500	kmer_match_length_7500_coverage_90pct	homolog_3
```

## R Interface

The homologs functionality is also available through R:

```r
library(alntools)

# Find homologous regions
results <- homologs_search(
  fasta_file = "assembly.fasta",
  query_contig = "ctg001",
  query_start = 1000,
  query_end = 11000,
  k = 21,
  num_kmers = 10,
  threshold = 80.0
)

# Results is a data.frame with additional columns:
# - length: Length of the homologous region
# - kmer_count: Number of matching kmers found
# - coverage_percent: Percentage of query kmers found
```
