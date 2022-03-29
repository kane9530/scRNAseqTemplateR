# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
#' auto_detect_tissue_type: automatically detect a tissue type of the dataset
#'
#' @param  path_to_db_file xlsx: DB file with cell types
#' @param scRNAseqData Matrix: input scRNA-seq matrix (rownames - genes, column names - cells),
#' @param seurat_object Seurat Object: Seurat object from which the scRNAseqdata is derived from
#' @param scaled Boolean: indicates whether the matrix is scaled. TRUE by default.
#' @importFrom openxlsx read.xlsx
#' @importFrom magrittr %>%
#' @importFrom grDevices rgb
#' @importFrom graphics barplot
#' @importFrom stats na.omit
#' @importFrom utils head
#' @import dplyr
#' @export

auto_detect_tissue_type <- function(path_to_db_file, scRNAseqData, seurat_object, scaled){

  # get all tissue types in DB
  db_read = read.xlsx(path_to_db_file); tissues_ = unique(db_read$tissueType); result_ = c()

  for(tissue in tissues_){ print(paste0("Checking...", tissue));

    # prepare gene sets
    gs_list = gene_sets_prepare(path_to_db_file, tissue);

    es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = scaled,
                            gs = gs_list$gs_positive, gs2 = gs_list$gs_negative);

    cL_resutls = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
    }))

    dt_out = cL_resutls %>% dplyr::group_by(.data$cluster) %>% dplyr::top_n(n = 1)

    # return mean score for tissue
    result_ = rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
  }

  # order by mean score
  result_ = result_[order(-result_$score),]

  # plot
  barplot(height=result_$score, names=result_$tissue, col=rgb(0.8,0.1,0.1,0.6),
          xlab="Tissue", ylab="Summary score",  main="The higher summary score, the more likely tissue type is")

  result_
}
