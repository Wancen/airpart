#' Plot allelic ratio as heatmap
#'
#' @param sce SingleCellExperiment
#' @param assay the assay to be plotted. Choices are \code{"ratio_pseudo"} which is the default,
#' \code{"ratio"}, \code{"counts"}.
#' @param genecluster an integer indicates which gene cluster heatmap want to be returned.
#' @param show_row_names show row names or not
#' @param order_by_group indicate whether order by group or order by cell types
#' @param ... Passsed on the other argument in
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return generates a heatmap
#'
#' @examples
#' set.seed(2021)
#' sce <- makeSimulatedData(p.vec = c(0.3, 0.5, 0.5, 0.3), ncl = 1)
#' sce <- preprocess(sce)
#' # display allelic ratio pattern in whole dataset
#' makeHeatmap(sce)
#'
#' sce <- geneCluster(sce, G = seq_len(4), plot = FALSE)
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' # display specific gene cluster partition result
#' makeHeatmap(sce_sub)
#' # display by cell type orders
#' makeHeatmap(sce_sub, order_by_group = FALSE)
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_block
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#'
#' @export
makeHeatmap <- function(sce, assay = c("ratio_pseudo", "ratio", "counts"), genecluster = NULL,
                        show_row_names = FALSE, order_by_group = TRUE, ...) {
  stopifnot("x" %in% names(colData(sce)))
  assay <- match.arg(assay, c("ratio_pseudo", "ratio", "counts"))
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  if (is.null(genecluster)) {
    m <- assays(sce)[[assay]]
  } else {
    m <- assays(sce[rowData(sce)$cluster == genecluster, ])[[assay]]
  }
  if (assay == "counts") {
    name <- "Total counts"
  } else {
    name <- "Allelic Ratio"
  }
  if (is.null(sce$part)) {
    split <- NULL
    ha <- ComplexHeatmap::HeatmapAnnotation(
      `cell type` = sce$x, border = FALSE,
      col = list(`cell type` = structure(col_vector[seq_len(nlevels(sce$x))],
        names = levels(sce$x)
      ))
    )
  } else {
    if (order_by_group) {
      split <- sce$part
      ha <- ComplexHeatmap::HeatmapAnnotation(
        group = anno_block(
          gp = gpar(fill = col_vector[seq_len(nlevels(split))]),
          labels_gp = gpar(col = "white", fontface = 4),
          labels = paste0("G", seq_len(nlevels(split)))
        ),
        `cell type` = sce$x, border = FALSE, group = sce$part,
        col = list(`cell type` = structure(col_vector[seq_len(nlevels(sce$x))],
          names = levels(sce$x)
        ))
      )
    } else {
      split0 <- rle(as.numeric(sce$part))
      y <- split0$lengths
      split <- rep(seq_along(y), y)
      ha <- ComplexHeatmap::HeatmapAnnotation(
        group = anno_block(
          gp = gpar(fill = col_vector[split0$values]),
          labels_gp = gpar(col = "white", fontface = 4),
          labels = paste0("G", split0$values)
        ),
        `cell type` = sce$x, border = FALSE,
        col = list(`cell type` = structure(col_vector[seq_len(nlevels(sce$x))],
          names = levels(sce$x)
        ))
      )
    }
  }
  ComplexHeatmap::Heatmap(m,
    column_split = split, column_title = NULL, top_annotation = ha,
    column_order = order(seq_len(ncol(m))), cluster_columns = FALSE,
    cluster_rows = FALSE, name = name,
    show_column_names = FALSE, show_row_names = show_row_names
  )
}
