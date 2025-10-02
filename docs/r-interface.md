# R Interface

`alntools` provides an R interface for constructing, loading, and querying alignment stores. The R interface offers the same functionality as the command-line tools but integrated into R workflows.

## Loading the R Interface

**CPU-only version:**
```R
library(Rcpp)
sourceCpp("cpp/aln_R.cpp")
```

## Available Functions

### Construction and Loading

```R
# Construct alignment store from PAF file
aln <- aln_construct(paf_file, max_reads = 0)

# Save alignment store to file
aln_save(aln, aln_file)

# Load existing alignment store
aln <- aln_load(aln_file)
```

**Parameters:**
- `paf_file`: Path to input PAF file
- `max_reads`: Maximum number of reads to process (0 = all)
- `aln_file`: Path to ALN file

### Query Functions

#### Bin Query

```R
bin_results <- aln_query_bin(aln, intervals, binsize, 
                            seg_threshold = 0.2, 
                            non_ref_threshold = 0.9, 
                            num_threads = 0, 
                            clip_mode_str = "all", 
                            clip_margin = 10, 
                            min_mutations_percent = 0.0, 
                            max_mutations_percent = 10.0,
                            min_alignment_length = 0, 
                            max_alignment_length = 0, 
                            min_indel_length = 3)
```

#### Pileup Query

```R
pileup_results <- aln_query_pileup(aln, intervals, report_mode, 
                                  clip_mode_str = "all", 
                                  clip_margin = 10, 
                                  min_mutations_percent = 0.0, 
                                  max_mutations_percent = 10.0,
                                  min_alignment_length = 0, 
                                  max_alignment_length = 0, 
                                  min_indel_length = 3)
```

**Parameters:**
- `report_mode`: "all", "covered", or "mutated"

#### Consensus Query

```R
consensus_results <- aln_query_consensus(aln, intervals, 
                                        consensus_threshold = 0.9, 
                                        min_consensus_coverage = 5, 
                                        num_threads = 0, 
                                        clip_mode_str = "all", 
                                        clip_margin = 10, 
                                        min_mutations_percent = 0.0, 
                                        max_mutations_percent = 10.0,
                                        min_alignment_length = 0, 
                                        max_alignment_length = 0, 
                                        min_indel_length = 3)
```

#### Full Query

```R
full_results <- aln_query_full(aln, intervals, height_style, 
                              max_alignments = 0, 
                              clip_mode_str = "all", 
                              clip_margin = 10,
                              min_mutations_percent = 0.0, 
                              max_mutations_percent = 10.0,
                              min_alignment_length = 0, 
                              max_alignment_length = 0, 
                              max_margin = 10, 
                              chunk_type, 
                              min_indel_length = 3)
```

**Parameters:**
- `height_style`: "by_coord_left", "by_coord_right", or "by_mutations"
- `chunk_type`: "read", "alignment", "break_on_overlap", or "break_on_gap"

**Returns:** List with `$alignments`, `$mutations`, `$reads`, and `$chunks` dataframes

#### Multi-library Variant Calling

```R
# Create a named list of alignment stores
store_list <- list("lib1" = aln1, "lib2" = aln2, "lib3" = aln3)

variants_results <- aln_query_variants(store_list, intervals, 
                                      min_variants_variant_support = 5, 
                                      min_variants_library_support = 2, 
                                      min_variants_coverage_support = 20, 
                                      clip_mode_str = "all", 
                                      clip_margin = 10,
                                      min_mutations_percent = 0.0, 
                                      max_mutations_percent = 10.0,
                                      min_alignment_length = 0, 
                                      max_alignment_length = 0, 
                                      min_indel_length = 3)
```

**Returns:** List with:
- `$variants`: Dataframe with variant information
- `$support`: Matrix with variant IDs as rows and library IDs as columns (read support)
- `$coverage`: Matrix with variant IDs as rows and library IDs as columns (coverage)
- `$library_ids`: Vector of library IDs

### Break Detection

```R
breaks_results <- aln_find_breaks(aln, window_size = 1000, p_threshold = 0.05, min_reads = 1)
```

**Returns:** Dataframe with columns: `contig`, `position`, `orientation`, `t`, `e`, `enrichment`, `pval`, `qval` (sorted by contig, then position)

## Common Parameters

**Alignment Filtering:**
- `clip_mode_str`: "all", "complete", "allow_one_side_clip", "only_one_side_clipped", "only_two_side_clipped", "only_clipped", "local_align"
- `clip_margin`: Margin in bases for clipping detection
- `min_mutations_percent`: Minimum mutation percentage
- `max_mutations_percent`: Maximum mutation percentage  
- `min_alignment_length`: Minimum alignment length
- `max_alignment_length`: Maximum alignment length (0 = no limit)
- `min_indel_length`: Minimum indel length for mutation calculations

**Threading:**
- `num_threads`: Number of threads (0 = auto-detect)

## Example R Script

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

# Run queries with filtering
bin_results <- aln_query_bin(aln, intervals, binsize, 
                            seg_threshold = 0.2, 
                            non_ref_threshold = 0.9, 
                            num_threads = 0,
                            min_alignment_length = 1000,  # Filter alignments < 1kb
                            min_indel_length = 5)         # Filter short indels < 5bp

pileup_results <- aln_query_pileup(aln, intervals, "covered", 
                                  max_alignment_length = 50000,  # Filter alignments > 50kb
                                  min_indel_length = 3)          # Filter short indels < 3bp

consensus_results <- aln_query_consensus(aln, intervals, 
                                        consensus_threshold = 0.8, 
                                        min_consensus_coverage = 10, 
                                        num_threads = 4,
                                        min_indel_length = 3)      # High-frequency variants ≥80% with ≥10x coverage

full_results <- aln_query_full(aln, intervals, "by_mutations", 
                              clip_mode_str = "complete",        # Only complete alignments
                              min_indel_length = 3)              # Filter short indels

breaks_results <- aln_find_breaks(aln, window_size = 1000, p_threshold = 0.05, min_reads = 3)

# Save results
write.table(bin_results, file = paste0(output_prefix, "_bins.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pileup_results, file = paste0(output_prefix, "_pileup.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(consensus_results, file = paste0(output_prefix, "_consensus.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_results$alignments, file = paste0(output_prefix, "_alignments.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_results$mutations, file = paste0(output_prefix, "_mutations.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_results$reads, file = paste0(output_prefix, "_reads.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(full_results$chunks, file = paste0(output_prefix, "_chunks.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(breaks_results, file = paste0(output_prefix, "_breaks.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
```

## Differences from CLI

The R interface provides the same core functionality as the command-line interface with these key differences:

1. **Rearrangement Analysis**: The R interface does not currently include rearrangement detection functionality. Use the CLI `rearrange` command for structural variation analysis.

2. **Verification**: Rearrangement verification is available in the R interface but is primarily used internally. The CLI `verify` command provides more comprehensive verification capabilities.

3. **Data Integration**: The R interface allows for seamless integration with R data analysis workflows, making it easier to combine alntools results with other analyses.

4. **Memory Management**: Alignment stores are kept in memory in R, which can be more efficient for multiple queries but may require more RAM for large datasets.

## Notes

- The R interface requires compilation of C++ code, which may take a few minutes on first use
- All output formats match those described in the [File Formats](file-formats.md) documentation
- Use the R interface when you need to integrate alignment analysis into larger R workflows
- For large-scale batch processing, the command-line interface may be more efficient
- The R interface supports the same filtering and parameter options as the CLI commands
