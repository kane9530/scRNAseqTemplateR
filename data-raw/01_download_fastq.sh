#!/bin/sh
# Download data from 10x genomics
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar -P ../inst/extdata

# Unzip data to recover fastq.gz files
cd ./inst/extdata
tar -xvf pbmc_1k_v3_fastqs.tar
ls 

# Remove .tar file
rm -rf pbmc*.tar
