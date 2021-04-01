#' Plot allelic ratio as heatmap
#'
#' @param sce SingleCellExperiment
#' @param show_row_names show row names or not
#' @param order_by_group indicate whether order by group or order by cell types
#' @param ... Passsed on the other argument in \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return generates a heatmap
#' 
#' @examples
#' sce <- makeSimulatedData(p.vec=c(0.3,0.5,0.3,0.3),ncl=1)
#' sce <- preprocess(sce)
#' # display allelic ratio pattern in whole dataset
#' makeRatioHeatmap(sce)
#'
#' sce <- geneCluster(sce, G = seq_len(4))
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' # display specific gene cluster partition result
#' makeRatioHeatmap(sce_sub)
#' # display by cell type orders
#' makeRatioHeatmap(sce_sub,order_by_group=FALSE)
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_block
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
makeRatioHeatmap <- function(sce, show_row_names = FALSE,
                             order_by_group = TRUE, ...) {
  m <- assays(sce)[["ratio_pseudo"]]
  if (is.null(sce$part)) {
    ha <- HeatmapAnnotation(
      `cell type` = sce$x, border = FALSE,
      col = list(`cell type` = structure(brewer.pal(nlevels(sce$x), "Set3"),
        names = levels(sce$x)
      ))
    )
    Heatmap(m,
      name = "Allelic Ratio", cluster_columns = FALSE, cluster_rows = FALSE,
      show_column_names = FALSE, show_row_names = show_row_names, top_annotation = ha
    )
  } else {
    if(order_by_group){
      split <- sce$part
      ha <- HeatmapAnnotation(
        group = anno_block(
          gp = gpar(fill = brewer.pal(9, "Pastel1")[seq_len(nlevels(split))]),
          labels_gp = gpar(col = "white",fontface=4),
          labels = paste0("group", seq_len(nlevels(split)))
        ),
        `cell type` = sce$x, border = FALSE,
        col = list(`cell type` = structure(brewer.pal(nlevels(sce$x), "Set3"), names = levels(sce$x)))
      )
      Heatmap(m,
              name = "Allelic Ratio", column_split = split,
              column_order = order(seq_len(ncol(m))),column_title = paste0("group", seq_len(nlevels(split))),
              cluster_columns = FALSE, cluster_rows = FALSE, show_column_names = FALSE,
              show_row_names = show_row_names, top_annotation = ha
      )
    }else{
      ha <- HeatmapAnnotation(
        group = sce$part, `cell type` = sce$x, border = FALSE,
        col = list(`cell type` = structure(brewer.pal(nlevels(sce$x), "Set3"),
                                           names = levels(sce$x)),
                   group = structure(brewer.pal(9, "Pastel1")[seq_len(nlevels(sce$part))],
                                     names = levels(sce$part)))
      )
      Heatmap(m,
              name = "Allelic Ratio", cluster_columns = FALSE, cluster_rows = FALSE,
              show_column_names = FALSE, show_row_names = show_row_names, top_annotation = ha)
    }

  }
}
