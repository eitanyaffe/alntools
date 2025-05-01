# load required libraries
library(ggplot2)

# plot alignments along contig, with their height, colored by mutation count
plot_alignments <- function(df, ofn) {
  # create plot
  p <- ggplot(df, aes(
    x = contig_start, xend = contig_end,
    y = height, yend = height,
    color = mutation_count
  )) +
    geom_segment(linewidth = 2) +
    scale_color_gradient(low = "#b6b6e5", high = "#ffbb00") +
    labs(x = "Contig Position", y = "Read Name", title = "Alignments along Contig") +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    )

  # save to file if specified
  if (!is.null(ofn)) {
    ggsave(ofn, p, width = 10, height = 8)
  }

  return(p)
}

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# check if we have the right number of arguments
if (length(args) < 1) {
  cat("usage: Rscript plot.r <alignment_file>\n")
  quit(status = 1)
}

# parse arguments
alignment_file <- args[1]
ofn <- args[2]

# read alignment data
df <- read.table(alignment_file, header = TRUE)

# generate and display/save plot
plot_alignments(df, ofn)
