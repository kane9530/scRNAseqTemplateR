#' Barcode rank / knee plot for filtering empty droplets
#'
#' Visualizes the inflection point to filter empty droplets.
#'
#' @param stats DataFrame: Output from `DropletUtil::barcodeRanks`
#' @param raw_cells Character vector: names of all unfiltered cells from raw_matrix
#' @param filt_cells  Character vector: names of filtered cells
#' @param save String: Output directory
#' @importFrom ggplot2 ggplot
#' @export

bc_rank_plot <- function(stats, raw_cells, filt_cells, save){
  Rx <- Tx <- cell <- NULL
  cells <- raw_cells %in% filt_cells
  keep <- !duplicated(stats$total)
  plot_df <- data.frame(Rx=stats$rank, Tx=stats$total, cell=cells)
  plot_df <- plot_df[keep, ]
  print({
    ggplot(subset(plot_df, plot_df$Tx >0), aes(x=Rx, y=Tx, col=cell, alpha=cell))+
      geom_point(size=3) +
      geom_hline(yintercept = stats@metadata$knee, lty = 2, col = '#0972D5', size=1.5) +
      annotate("text", x=max(plot_df$Tx), y=stats@metadata$knee+10000, label="Knee", color = "#0972D5", size=5.5) +
      geom_hline(yintercept = stats@metadata$inflection, lty = 2, col = '#09CFD5', size = 1.5) +
      annotate("text", x=max(plot_df$Tx), y=stats@metadata$inflection+500, label="Inflection", color = "#09CFD5", size=5.5) +
      scale_x_log10(labels = plain, breaks = scales::trans_breaks("log10", function(x) round(10^x, 0))) +
      scale_y_log10(breaks = scales::trans_breaks('log10', function(x) floor(10^x)), labels = plain)+
      scale_color_manual(values = c('#8595A8', '#6406B6'), name = NULL, labels = c("Background", "Cells")) +
      scale_alpha_manual(values = c(0.5,1)) +
      labs(x = 'Barcodes', y = 'UMI counts', title = 'Barcode Rank Plot') +
      guides(alpha = "none", colour = guide_legend(reverse = TRUE, override.aes=list(size = 5))) +
      theme_linedraw() +
      theme(plot.title = element_text(size=20, hjust=0.5, face = 'bold'),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 15),
            legend.text = element_text(size=19),
            legend.background = element_rect(fill = 'transparent'),
            legend.position = c(0.15,0.15))
  })
}
