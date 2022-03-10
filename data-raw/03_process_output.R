# Run this after running 01_download_fastq.sh and 02_kallisto-bustools.sh sequentially
# This file imports the kallisto-output data into R

# Generate raw_matrix.rda
raw_matrix <- read_count_output("../inst/extdata/output/counts_unfiltered", name = "cells_x_genes")

# Generate tr2g.rda
tr2g <- readr::read_tsv("../inst/extdata/t2g.txt", col_names = c("transcript", "gene", "gene_name"))

