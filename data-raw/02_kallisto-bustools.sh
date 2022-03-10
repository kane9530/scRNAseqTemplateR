#!/bin/sh
# Download a pre-built kallisto transcriptome index (We are not building this locally)
kb ref -d human -i ../inst/extdata/index.idx -g ../inst/extdata/t2g.txt -f1 ../inst/extdata/transcriptome.fasta

# Pseudoalignment and transcript quantification via kallisto-bustools
# Run command in directory containing the index.idx and t2g.txt files, as well as the pbmc_1k_v3_fastqs directory
kb count -i index.idx -g t2g.txt -x 10xv3 -o output --filter bustools -t 2 ../inst/extdata/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz ../inst/extdata/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz ../inst/extdata/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz ../inst/extdata/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz
