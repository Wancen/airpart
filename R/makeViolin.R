#' Plot allelic ratio as heatmap
#'
#' @param sce SingleCellExperiment
#'
#' @return a ggplot2 object
#'
#' @examples
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' makeViolin(sce_sub)
#' @importFrom ggplot2 ggplot aes geom_violin scale_fill_brewer theme labs
#'
#' @export
makeViolin <- function(sce) {
  cl_ratio <- as.vector(unlist(assays(sce)[["ratio"]]))
  cl_total <- as.vector(unlist(counts(sce)))
  dat <- data.frame(
    ratio = cl_ratio,
    x = factor(rep(sce$x, each = length(sce))),
    cts = cl_total,
    part = factor(rep(sce$part, each = length(sce)))
  )
  dat <- dat[!is.nan(dat$ratio), ]
  p <- ggplot(dat, aes(x = .data$x, y = .data$ratio, fill = .data$part)) +
    geom_violin(alpha = 0.5) +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "BuPu") +
    labs(x = "cell type", y = "allelic ratio")
  p
}
