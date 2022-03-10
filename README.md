
# scRNAseqTemplateR

scRNAseqTemplateR was created as a template for the analysis of single-cell RNA sequencing (scRNA-seq) data in the `Genomics and Data Analytics Core (GeDaC)`, Cancer Science Institute of Singapore. 

The package organises the analysis of the small 10x Chromium `PBMC_V3_1K` dataset into a reproducible and portable fashion. Other scRNA-seq analysis pipelines can be adapted from this template. The R package facilitates bioinformaticians and biologists to share and reuse the code, and to tweak and extend the analysis as desired.

## Installation

You can install the development version of scRNAseqTemplateR like so:

``` r
devtools::install_github("kane9530/scRNAseqTemplateR")

```

## Deliverables

1. A summary of the methods employed can be found in `results/methods.pdf`. This is suitable for a direct copy-paste to the `Methods` section of a publication.

2. Results files (e.g. `.csv`) and plots (`.png`) can be located in the `results/` folder.

3. A report of the analysis steps taken can be found in static `.pdf` or interactive `.html` formats in `vignettes/`.

## Directories

- DESCRIPTION: gives an overview of the project and its dependencies.

- data/ contains the .rda tidy data files

- data-raw/ contains the `. sh` and `.R` scripts to download the fastq files, pseudoalign reads with kallisto-bustools and generate the `raw_matrix.rda` and `t2g.rda` files in /data

- man/ contains the documentation for the data and functions

- vignettes/ contains the analysis report  as a package vignette

- R/: contains R scripts with functions used throughout the package

- tests/: contains development-time tests for our functions

- results/ contains all the outputs from the vignette, as well as a summary of the methods used in the analysis

## Examples

``` r
# Browse the vignette here:

# Load the data

# Visualising plot from seurat object



```

