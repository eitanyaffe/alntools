# Test script for querying alignments by read IDs

# Load the Rcpp package
library(Rcpp)

# Compile and load the C++ code
cdir <- ".Rcpp_dir"
sourceCpp("cpp/aln_R.cpp", verbose = TRUE, cacheDir = cdir)

# Input alignment file path
# Modify this to point to your actual alignment file
aln_file <- "examples/test.aln"

# Check if the file exists
if (!file.exists(aln_file)) {
  stop(paste("Alignment file not found:", aln_file))
}

# Load the alignment
cat("Loading alignment from:", aln_file, "\n")
aln <- aln_load(aln_file)

# Define a list of read IDs to query
# Replace these with actual read IDs from your data
read_ids <- c("read1", "read2", "read3")
cat("Querying for read IDs:", paste(read_ids, collapse = ", "), "\n")

# Perform the query
result <- aln_alignments_from_read_ids(aln, read_ids)

# Print query results
cat("Found", nrow(result), "alignments for specified read IDs\n")
if (nrow(result) > 0) {
  # Print the first few rows
  cat("First few results:\n")
  print(head(result))

  # Count alignments per read
  read_counts <- table(result$read_id)
  cat("\nAlignments per read:\n")
  print(read_counts)
}

# Save to file if desired
output_file <- "read_query_results.tsv"
write.table(result, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Results saved to:", output_file, "\n")
