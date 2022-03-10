#' Plot pct genes
#'
#' Plot genes with the highest proportions cells. Adapted from kallistobustools tutorial
#'
#' @param mat Matrix/dGCmatrix: Counts matrix
#' @param tr2g Dataframe/tibble: Tr2g output from kb workflow, which maps ensembl transcript IDs to gene names and gene IDs.
#' @param top_n Integer: top n differentially expressed genes
#' @param symbol string: rownames/gene names. Can be either "ensembl" or "symbol".
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom forcats fct_reorder
#' @export

plot_pct_genes <- function(mat, tr2g, top_n = 20, symbol = "ensembl") {
  gene <- gene_name <- value <- NULL
  # Total gene count, summed over all cells (barcodes)
  pct_tx <- Matrix::rowSums(mat)
  # List of gene names in descending order of counts
  gs <- rownames(mat)[order(-pct_tx)]
  # Df of Cells(rows) x genes(cols) of top_n most highly expressed genes
  df <- as.data.frame(t(as.matrix(mat[gs[1:top_n],])))
  # Normalise all gene counts by total number of counts in cells
  # Obtain a longer #gene x 2 df where all genes collapsed into 1 'gene' column
  df <- df %>%
    dplyr::mutate_all(function(x) x/Matrix::colSums(mat)) %>%
    tidyr::pivot_longer(everything(), names_to = "gene")

  if (symbol == "ensembl") {
    df <- dplyr::left_join(df, tr2g, by = "gene")
  }
  else if (symbol == "symbol"){
    df <- dplyr::rename(df, gene_name = gene)
  }
  df %>%
    dplyr::mutate(gene = forcats::fct_reorder(gene_name, value, .fun = stats::median)) %>%
    ggplot(ggplot2::aes(gene, value)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(x = "", y = "Proportion of total counts") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal()
}
