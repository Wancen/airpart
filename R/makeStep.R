#' Plot group partition and posterior allelic ratio estimates by step
#'
#' @param sce SingleCellExperiment
#' @param xlab the x axis name.
#'
#' @return a ggplot2 object.
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' sce_sub <- allelicRatio(sce_sub)
#' makeStep(sce_sub)
#' @export
makeStep <- function(sce, xlab = "cell type") {
  ar <- rowData(sce)[, c(grep("ar_", colnames(rowData(sce)), value = TRUE))] %>%
    `colnames<-`(levels(sce$x))
  dat <- data.frame(
    ratio = as.vector(unlist(ar)),
    x = factor(rep(unique(sce$x), each = length(sce))),
    part = factor(rep(metadata(sce)$partition$part, each = length(sce))),
    feat = factor(rep(row.names(sce), nlevels(sce$x)), levels = row.names(sce))
  )
  p <- ggplot2::ggplot(dat, aes(x = .data$x, y = .data$ratio)) +
    ggplot2::geom_point(aes(color = .data$part)) +
    ggplot2::geom_step(group = 1) +
    ggplot2::geom_hline(yintercept = 0.5, colour = "gray40", linetype = "dashed") +
    ggplot2::facet_wrap(~feat) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = "transparent")
    ) +
    labs(x = xlab, y = "allelic ratio")
  p
}
