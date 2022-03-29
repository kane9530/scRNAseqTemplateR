## ----chunk_options, include = FALSE-------------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.width = 7, fig.height = 5)

## ---- load_packages,  include=FALSE-------------------------------------------
library(scRNAseqTemplateR)
library(dplyr)

## ----loadData, include=FALSE--------------------------------------------------
# Raw count_matrix - Output of Kallisto-bustools
load(file = "../data/raw_matrix.rda")

# Transcript to genes dataframe - Output of Kallisto-bustools
load(file= "../data/tr2g.rda")

# Defining directories

## Kallisto output directory
## This directory contains the output files from the kallisto-bustools pseudoalignment steps, and are stored in inst/extdata
kallisto_path <- '../inst/extdata/output'

## Saving all results
results_path = '../results'

## Saving the filtered matrix in 10x format 
out_path = '../inst/extdata/output_filt/counts_filtered'


## ----UMIcounts_barcodes-------------------------------------------------------
# Dimensions of matrix
dim(raw_matrix)

# Examining distribution of UMI counts over barcodes
tot_counts <- Matrix::colSums(raw_matrix)
summary(tot_counts)

## ----Plot_genes_proportion----------------------------------------------------
plot_pct_genes(raw_matrix, tr2g)

## ----removeEmpty, results=FALSE-----------------------------------------------
# Use DropletUtils package to get probability that each barcode is a cell
out <- DropletUtils::emptyDrops(raw_matrix)

# Set threshold probability for calling a cell
keep <- out$FDR <= 0.05 

# Use threshold to remove empty drops
keep[is.na(keep)] <- FALSE
filt_matrix <- raw_matrix[,keep] 

# Dimensions of filtered matrix: 60623 x 1222
head(filt_matrix)

# Proportion of cells kept: Approx 2%
dim(filt_matrix)[2]/dim(raw_matrix)[1]*100

# Keep distinct rows of tr2g

tr2g <- dplyr::distinct(tr2g[, c("gene", "gene_name")])


## ----write10xcounts_filteredMat, eval=FALSE-----------------------------------
#  
#  DropletUtils::write10xCounts(out_path,
#                 x = filt_matrix,
#                 gene.id = unlist(tr2g[,1]),
#                 gene.symbol = unlist(tr2g[,2]),
#                 overwrite=T)

## ----calculate-stats, results=FALSE-------------------------------------------

# Load filtered mtx
#filt_mtx <- Matrix::readMM(paste0(out_path,'/matrix.mtx')) 
load(file="../data/filt_matrix.rda")

# Load run info from JSON files produced by Kb
#kb_stats <- c(fromJSON(file = '../data/output/inspect.json'), 
#              fromJSON(file = '../data/output/run_info.json')) 
load(file = "../data/kb_stats.rda")

# Determine 10x Chromium library prep chemistry version
tech <- grep('10[X/x](.*)', strsplit(kb_stats$call, '\\s')[[1]], value=T) 

# Summarise stats in 'pretty' way
seq_stats <- data.frame(stat = c('Sequencing technology', 
                                 'Number of reads processed', 
                                 '% reads pseudoaligned', 
                                 # get sequencing/alignment stats 
                                 '% reads on whitelist'), 
                        value = prettyNum(c(tech, 
                                kb_stats$n_processed,
                                kb_stats$p_pseudoaligned,
                                round(kb_stats$percentageReadsOnWhitelist,2)), 
                                big.mark = ','))

# Calculate cell stats and save to df
# Total percentage of counts in filtered vs raw mtx
p_cnts_in_cells <- round((sum(filt_matrix)/sum(raw_matrix))*100, 2) 
# Median number of gene counts in cells in filtered matrix
med_cnts_cell <- median(colSums(filt_matrix))
# Median number of genes expressed (>=1 counts) in filtered matrix
med_genes_cell <- median(apply(filt_matrix, 2, function(x) sum(x >= 1)))
# Total number of genes detected
tot_genes_detected <- sum(rowSums(filt_matrix)>=1)

cell_stats <- data.frame(stat = c('Estimated number of cells', 
          '% counts in cells', 
          'Median counts per cell', 
          'Median genes per cell', 
          'Total genes detected'), 
          value = prettyNum(c(ncol(filt_matrix), p_cnts_in_cells, med_cnts_cell,
          med_genes_cell, tot_genes_detected), big.mark = ','))

# Get rank stats and save
stats <- DropletUtils::barcodeRanks(raw_matrix)

#write.table(stats, file=paste0(out_path, "/barcode_rank.txt"))


## ----knee-plot----------------------------------------------------------------
# load raw cells
#raw_cells <- read.table(paste0(kallisto_path,'/counts_unfiltered/cells_x_genes.barcodes.txt'),header=FALSE)[,1] 
load(file="../data/raw_cells.rda")

# load filtered cells
#filt_cells <- read.csv(paste0(out_path, "/barcodes.tsv"), header = F, sep ='\t')[,1] 
load(file="../data/filt_cells.rda")

# create barcode rank plot and save as .png
bc_rank_plot(stats = stats, raw_cells = raw_cells, filt_cells = filt_cells, 
             save ="/plots/barcode_rank.png")


## ----createSeuratObject, eval = FALSE-----------------------------------------
#  # Create expression matrix with Read10X function from Seurat
#  # out_path contains the genes.tsv, barcodes.tsv and matrix.mtx files
#  # gene.column = 2 means we use the gene_symbols (not the ensembl IDs) for gene names
#  
#  expression_matrix <- Read10X(
#    out_path,
#    gene.column = 2,
#    cell.column = 1,
#    unique.features = TRUE
#  )
#  
#  # Create Seurat object
#  pbmc.1k.seurat <- CreateSeuratObject(counts = expression_matrix,
#                                       project = "pbmc1k",
#                                       min.cells = 3,
#                                       min.features = 200)
#  

## ----load_from_data, include=FALSE--------------------------------------------

load(file = "../data/pbmc.1k.seurat.rda")


## ---- QC_pert_Mito------------------------------------------------------------
# Stash the percent.mt stats in seurat object metadata
# NOTE: Change 'MT' to 'mt' for mouse
pbmc.1k.seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = pbmc.1k.seurat, pattern = "^MT") 

# Violin plots of the 3 QC metrics
Seurat::VlnPlot(pbmc.1k.seurat, c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        pt.size = 0.1, ncol=3)

## ---- QC_scatter_genes_mito---------------------------------------------------

p1 <- ggplot2::ggplot(pbmc.1k.seurat@meta.data, ggplot2::aes(nCount_RNA, nFeature_RNA)) +
  ggplot2::geom_point(alpha = 0.7, size = 0.5) +
  ggplot2::labs(x = "Total UMI counts per cell", y = "Number of genes detected")+
  ggplot2::theme_minimal()

p2 <- ggplot2::ggplot(pbmc.1k.seurat@meta.data, ggplot2::aes(nCount_RNA, percent.mt)) +
  ggplot2::geom_point(alpha = 0.7, size = 0.5) +
  ggplot2::labs(x = "Total UMI counts per cell", y = "Percentage of mitochondrial genes")+
  ggplot2::theme_minimal()

# 1x2 layout
gridExtra::grid.arrange(p1, p2, ncol=2)


## ----subset_seurat------------------------------------------------------------
pbmc.1k.seurat <-subset(pbmc.1k.seurat, subset = 
                           nCount_RNA < 20000 & 
                           nCount_RNA > 1000 & 
                           nFeature_RNA > 1000 & 
                           percent.mt < 20)

## ----QC_top_n_counts----------------------------------------------------------
verbose = TRUE
pbmc.1k.seurat.2 <- pbmc.1k.seurat %>%
  Seurat::NormalizeData(verbose = verbose) 

plot_pct_genes(Seurat::GetAssayData(pbmc.1k.seurat.2, slot = "counts"), 
               tr2g, symbol = "symbol")

## ----top2k_var_genes----------------------------------------------------------
pbmc.1k.seurat.3 <- pbmc.1k.seurat.2 %>%
  Seurat::FindVariableFeatures(verbose = verbose) 

# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(pbmc.1k.seurat.3), 10)

# plot variable features with and without labels
plot1 <- Seurat::VariableFeaturePlot(pbmc.1k.seurat.3, log = FALSE)
Seurat::LabelPoints(plot = plot1, points = top10, repel = FALSE)


## ----scale_pca, fig.width = 7, fig.height = 5---------------------------------
pbmc.1k.seurat.4 <- pbmc.1k.seurat.3 %>%
  Seurat::ScaleData(verbose = verbose) %>%
  Seurat::RunPCA(npcs = 40, verbose = verbose) 

Seurat::VizDimLoadings(pbmc.1k.seurat.4, dims= 1:2, reduction = "pca")

## ----visualise_gene-----------------------------------------------------------
Seurat::FeaturePlot(pbmc.1k.seurat.4, reduction = "pca", feature = "CST3")

## ----elbowplot----------------------------------------------------------------
Seurat::ElbowPlot(pbmc.1k.seurat.4)


## ----louvain_clustering-------------------------------------------------------
pbmc.1k.seurat.5<- pbmc.1k.seurat.4 %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:10, verbose = verbose) %>% 
    Seurat::FindClusters(resolution = 0.5, verbose = verbose)

## ----umap_embedding-----------------------------------------------------------
pbmc.1k.seurat.6<- pbmc.1k.seurat.5 %>%
  Seurat::RunUMAP(reduction = "pca", dims = 1:10, verbose = verbose)

Seurat::DimPlot(pbmc.1k.seurat.6, reduction = "umap", split.by = "orig.ident", label = TRUE)


## ----umap_gene_expression-----------------------------------------------------
Seurat::FeaturePlot(pbmc.1k.seurat.6, reduction = "umap", features = c("CST3", "NKG7", "PPBP"),
ncol = 3)

## ----findmarkers--------------------------------------------------------------
top_n = 20
cluster1.markers <- Seurat::FindMarkers(pbmc.1k.seurat.6, ident.1 = 1, min.pct = 0.25)
cluster1.markers$pct.diff <- cluster1.markers$pct.1 - cluster1.markers$pct.2
cluster1.markers.df <- as_tibble(cluster1.markers, rownames = "geneID")

# Export DEGs for each cluster (ranked by avg_logFC > 0.5)
myTopHits_cluster1 <- cluster1.markers.df %>% arrange(desc(avg_log2FC))
myTopHits_cluster1 <- dplyr::slice(myTopHits_cluster1, 1:top_n)
#myTopHits_cluster1
# Interactive table 
DT::datatable(myTopHits_cluster1, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: Cluster 1 genes',
          options = list(keys = TRUE, searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>%
  DT::formatRound(columns=c(2:ncol(myTopHits_cluster1)), digits=2)

## ----findallmarkers, fig.width = 9, fig.height = 10---------------------------
pbmc.1k.markers <- Seurat::FindAllMarkers(pbmc.1k.seurat.6, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)

# Top 10 marker genes for each cluster and plot as a heatmap
top3 <- pbmc.1k.markers %>% 
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 3, wt = avg_log2FC)

Seurat::DoHeatmap(pbmc.1k.seurat.6, features = top3$gene)

## ----export_markers, eval=FALSE-----------------------------------------------
#  filename = paste0(results_path,"/markers_all.csv")
#  write.csv(pbmc.1k.markers %>% group_by(cluster), file=filename)
#  

## ----featuremap_allmarkers, fig.width = 9, fig.height = 10--------------------

Seurat::FeaturePlot(pbmc.1k.seurat.6, features = top3$gene, ncol = 5)


## ----violinplots_allmarkers, fig.width = 7, fig.height = 10-------------------

Seurat::VlnPlot(pbmc.1k.seurat.6, features = top3$gene, ncol = 5)


## ----dotplot_manualAnnotate---------------------------------------------------

Seurat::DotPlot(pbmc.1k.seurat.6, assay = "RNA", features = top3$gene, scale.by = "size") +
  ggplot2::coord_flip()


## ----umap_with_ident----------------------------------------------------------

new.cluster.ids <- c("CD14+ Mono", "Memory CD4 T", "Naive CD4 T", "B1", "FCGR3A+ Mono", 
    "NK", "CD8+ T", "B2", "Platelet")
names(new.cluster.ids) <- levels(pbmc.1k.seurat.6)
pbmc.1k.seurat.6 <- Seurat::RenameIdents(pbmc.1k.seurat.6, new.cluster.ids)
Seurat::DimPlot(pbmc.1k.seurat.6, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4)


## ----singleR_heatmap----------------------------------------------------------

Monaco.data <- celldex::MonacoImmuneData(ensembl = FALSE) 

# Converting from seuratobject to singlecellexperiment container
pbmc.1k.sce <- Seurat::as.SingleCellExperiment(pbmc.1k.seurat.6)

# Leverage the Monaco dataset for prediction
# Can change the ref data to annotate with different datasets

predictions <- SingleR::SingleR(test=pbmc.1k.sce, assay.type.test=1, 
                       ref=Monaco.data, labels=Monaco.data$label.main)

SingleR::plotScoreHeatmap(predictions)

## ----singleR_predictions, fig.width = 7, fig.height = 10----------------------

#Add prediction labels back to singleCellExperiment object
pbmc.1k.sce[["SingleR.labels"]] <- predictions$labels

#Convert back to seurat object to utilise the seurat dimPlot function
pbmc.1k.seurat.7 <- Seurat::CreateSeuratObject(counts = pbmc.1k.sce@assays@data@listData[["counts"]],
                                       meta.data = as.data.frame(colData(pbmc.1k.sce)),
                                       reduction = pbmc.1k.sce@int_colData$reducedDims$UMAP)

pbmc.1k.seurat.7[["UMAP"]] <- Seurat::CreateDimReducObject(embeddings = pbmc.1k.sce@int_colData$reducedDims$UMAP, assay = Seurat::DefaultAssay(pbmc.1k.seurat.7))

#Plot UMAP
Seurat::DimPlot(pbmc.1k.seurat.7, reduction = "UMAP", 
        group.by = "SingleR.labels", 
        label = TRUE,
        repel = FALSE)+ 
  patchwork::plot_annotation(title = 'Labelling clusters with singleR labels')+
  ggplot2::theme(plot.title = ggplot2::element_blank())


## ----session info, include=TRUE-----------------------------------------------
sessionInfo()

