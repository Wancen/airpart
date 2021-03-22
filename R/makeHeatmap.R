#' Plot allelic ratio as heatmap
#'
#' @param sce SingleCellExperiment
#' @param show_row_names show row names or not
#' @param ... Passsed on the other argument in \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @examples
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' # display allelic ratio pattern in whole dataset
#' makeRatioHeatmap(sce)
#'
#' sce <- geneCluster(sce, G=1:4)
#' sce_sub <- wilcoxExt(sce,genecluster=1)
#' # display specific gene cluster partition result
#' makeRatioHeatmap(sce_sub)
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_block
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
makeRatioHeatmap <- function(sce, show_row_names = FALSE, ...) {
  m <- assays(sce)[["ratio_pseudo"]]
  if (is.null(sce$part)) {
    ha <- HeatmapAnnotation(`cell type` = sce$x, border = FALSE,
                            col = list(`cell type` = structure(brewer.pal(nlevels(sce$x), "Set3"),
                                                               names = levels(sce$x))))
    Heatmap(m, name = "Allelic Ratio", cluster_columns = FALSE, cluster_rows = FALSE,
            show_column_names = FALSE, show_row_names = show_row_names, top_annotation = ha)
  } else {
    split <- sce$part
    ha <- HeatmapAnnotation(
      group = anno_block(
        gp = gpar(fill = brewer.pal(9, "Pastel1")[1:nlevels(split)]),
        labels_gp = gpar(col = "white", fontsize = 10),
        labels = paste0("group", 1:nlevels(split))
      ),
      `cell type` = sce$x, border = FALSE,
      col = list(`cell type` = structure(brewer.pal(nlevels(sce$x), "Set3"), names = levels(sce$x)))
    )
    Heatmap(m, name = "Allelic Ratio", column_split = split,
            cluster_columns = FALSE, cluster_rows = FALSE, show_column_names = FALSE,
            show_row_names = show_row_names, top_annotation = ha, ...)
  }
}
