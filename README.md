# Extracting_markers using R-programming
# Install and load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

library(Biostrings)


setwd("/home/Desktop/markers/chr7")
getwd()

# Define the path to your FASTA file
fasta_file <- "chr7.fa"

# Read the FASTA file
chr7_seq <- readDNAStringSet(fasta_file)

# Extract the sequence
seq <- chr7_seq[[1]]

# Get the sequence length
seq_length <- length(seq)
cat("Length of the chromosome 7 sequence:", seq_length, "\n")

# Define marker size and step size
marker_size <- 2000
step_size <- 100000

# Create a list to store markers
markers <- list()

# Extract markers every 100,000 base pairs, with each marker being 2000 bases
for (start_pos in seq(1, seq_length - marker_size + 1, by = step_size)) {
  marker <- subseq(seq, start = start_pos, end = start_pos + marker_size - 1)
  markers[[length(markers) + 1]] <- list(
    start = start_pos,
    end = start_pos + marker_size - 1,
    sequence = marker
  )
}

cat("Number of markers extracted:", length(markers), "\n")

# Save markers in FASTA format with the desired header format
output_fasta_file <- "new_markers_chr7.fa"
fasta_lines <- character()

for (i in seq_along(markers)) {
  header <- paste0(">chr7:", markers[[i]]$start, "-", markers[[i]]$end, " marker", i)
  sequence <- as.character(markers[[i]]$sequence)
  fasta_lines <- c(fasta_lines, header, sequence)
}

writeLines(fasta_lines, output_fasta_file)
cat("Markers saved to:", output_fasta_file, "\n")
