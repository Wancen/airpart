#' Plot allelic ratio as heatmap
#'
#' @param se SummarizedExpeirment
#'
#' @importFrom pheatmap pheatmap
#' @export
plotRatioHeatmap <- function(se) {
  anno_df <- data.frame(x=factor(se$x), row.names=colnames(se))
  pheatmap(assays(se)[["ratio_pseudo"]],
           color=colorRampPalette(c("blue","white","red"))(101),
           breaks=0:100/100,
           annotation_col=anno_df,
           cluster_rows=FALSE, cluster_cols=FALSE,
           show_colnames=FALSE, show_rownames=FALSE)
}
