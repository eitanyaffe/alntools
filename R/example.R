# Load the Rcpp package
library(Rcpp)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument (the .aln file) is provided
if (length(args) < 1) {
  stop("Usage: Rscript example.R <path_to_paf_file> <path_to_intervals_file> <binsize> <ofn_prefix>", call. = FALSE)
}

paf_ifn <- args[1]
intervals_ifn <- args[2]
binsize <- as.numeric(args[3])
ofn_prefix <- args[4]

# print the arguments
cat(paste0("paf_ifn: ", paf_ifn, "\n"))
cat(paste0("intervals_ifn: ", intervals_ifn, "\n"))
cat(paste0("binsize: ", binsize, "\n"))
cat(paste0("ofn_prefix: ", ofn_prefix, "\n"))

# Compile and load the C++ code (can specify build directory)
cdir <- ".Rcpp_dir"
rebuild <- F
sourceCpp("cpp/aln_R.cpp", verbose = T, cacheDir = cdir, rebuild = rebuild)

################################################################################
# Construction example
################################################################################

construct <- function(paf_file, fn, max_reads = 0) {
  cat("constructing alignment from PAF example\n")
  # Create the alignment store
  aln <- aln_construct(paf_file, max_reads)

  # Save the alignment to file
  cat("saving alignment to", fn, "\n")
  aln_save(aln, fn)
}

################################################################################
# QueryBin
################################################################################

query_bins <- function() {
  cat("querying bin example\n")
  bin_results <- aln_query_bin(aln, intervals, binsize)
  ofn <- paste0(ofn_prefix, "_bins.tsv")
  cat(paste0("saving output to ", ofn, "\n"))
  write.table(bin_results, file = ofn, sep = "\t", row.names = F, quote = F)
}

################################################################################
# QueryPileup
################################################################################

query_pileup <- function() {
  cat("querying pileup example\n")
  pileup_results <- aln_query_pileup(aln, intervals, "covered")
  ofn <- paste0(ofn_prefix, "_pileup.tsv")
  cat(paste0("saving output to ", ofn, "\n"))
  write.table(pileup_results, file = ofn, sep = "\t", row.names = F, quote = F)
}

################################################################################
# QueryFull
################################################################################

query_full <- function() {
  cat("querying full example\n")
  full_results <- aln_query_full(aln, intervals, "by_mutations")
  ofn_alignments <- paste0(ofn_prefix, "_alignments.tsv")
  ofn_mutations <- paste0(ofn_prefix, "_mutations.tsv")
  cat(paste0("saving alignments to ", ofn_alignments, "\n"))
  cat(paste0("saving mutations to ", ofn_mutations, "\n"))
  write.table(full_results$alignments, file = ofn_alignments, sep = "\t", row.names = F, quote = F)
  write.table(full_results$mutations, file = ofn_mutations, sep = "\t", row.names = F, quote = F)
}

################################################################################
# main
################################################################################

tryCatch(
  {
    aln_fn <- paste0(ofn_prefix, ".aln")
    construct(paf_ifn, aln_fn)
    aln <- aln_load(aln_fn)
    intervals <- read.table(intervals_ifn, header = T)
    cat("number of input intervals: ", nrow(intervals), "\n")

    query_bins()
    query_pileup()
    query_full()
  },
  error = function(e) {
    cat("Error: ", e$message, "\n")
  }
)
