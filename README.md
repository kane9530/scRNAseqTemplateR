# scRNAseqTemplateR

The `scRNAseqTemplateR` R package was created to provide a portable and reproducible template for the analysis of single-cell RNA sequencing (scRNA-seq) data. The R package provides a starting template for bioinformaticians and biologists to share and reuse the code, and to adapt and extend the analysis as required.

To illustrate, we present the analysis of the 10x Chromium `PBMC_V3_1K` dataset. Other scRNA-seq analysis pipelines can be easily adapted from this template. 

## Installation

You can install the development version of scRNAseqTemplateR like so:

``` r
devtools::install_github("kane9530/scRNAseqTemplateR")

```

## Deliverables

Within the package, all output files are located in the`/results` and `/vignettes` directories. Output files can be classified into 2 categories:

### Quick summaries

1. A summary of the overall results can be found in `results/summary.pdf`. This consolidates the critical quality control assessments and downstream analyses into  quick, bite-sized, jargon-free reading.

2. A description the methods employed can be found in `results/methods.pdf`. This is suitable for a direct copy-paste to the `Methods` section of a publication.

### Detailed results and report 

1. Results files (e.g. `.csv`) and plots (`.png`) can be located in the `results/` folder.

2. A report of the analysis steps taken can be found in static `.pdf` or interactive `.html` formats in `vignettes/`.

## Directories

- `DESCRIPTION`: gives an overview of the project and its dependencies.

- `data/` contains the .rda tidy data files

- `data-raw/` contains the `. sh` and `.R` scripts to download the fastq files, pseudoalign reads with kallisto-bustools and generate the `raw_matrix.rda` and `t2g.rda` files in /data

- `man/` contains the documentation for the data and functions

- `vignettes/` contains the analysis report  as a package vignette

- `inst/doc` contains vignette files copied to the top-level directory

- `R/`: contains R scripts with functions used throughout the package

- `tests/`: contains development-time tests for our functions

- `results/` contains all the outputs from the vignette, as well as a summary of the methods used in the analysis

## Usage

``` r
library(scRNAseqTemplateR)

# Browse the vignette here:
browseVignettes("scRNAseqTemplateR")

# Load the PBMC data - Loadable data can be found in /data
data(raw_matrix)

# Checking for function documentation 
?plot_pct_genes

```

## Acknowledgements

The work was performed in the `Genomics and Data Analytics Core (GeDaC)`, Cancer Science Institute of Singapore. 
