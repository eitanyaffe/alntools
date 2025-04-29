# Load the Rcpp package
library(Rcpp)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument (the .aln file) is provided
if (length(args) < 1) {
  stop("Usage: Rscript example.R <path_to_aln_file> <path_to_intervals_file> <binsize> <ofn_prefix>", call. = FALSE)
}

aln_ifn <- args[1]
intervals_ifn <- args[2]
binsize <- as.numeric(args[3])
ofn_prefix <- args[4]

# print the arguments
cat(paste0("aln_ifn: ", aln_ifn, "\n"))
cat(paste0("intervals_ifn: ", intervals_ifn, "\n"))
cat(paste0("binsize: ", binsize, "\n"))
cat(paste0("ofn_prefix: ", ofn_prefix, "\n"))

# Compile and load the C++ code (can specify build directory)
cdir <- ".Rcpp_dir"
rebuild <- F
sourceCpp("cpp/aln_R.cpp", verbose = T, cacheDir = cdir, rebuild = rebuild)

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

tryCatch(
  {
    aln <- aln_load(aln_ifn)
    intervals <- read.table(intervals_ifn, header = T)
    cat("number of input intervals: ", nrow(intervals), "\n")

    query_bins()
  },
  error = function(e) {
    cat("Error: ", e$message, "\n")
  }
)
