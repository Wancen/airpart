#' Estimate overdispersion parameter of a betabinomial
#'
#' @param sce SingleCellExperiment with \code{ase.mat} and \code{counts}
#' @param ct the name of the cell type to be estimated. Default is the highest one.
#' @param pc pseudocount for calculating the smoothed ratio in the preprocess step.
#' @param genecluster which gene cluster dispersion parameter want to be estimated.
#' Default is the cluster with the most cells
#'
#' @examples
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' estDisp(sce)
#' @importFrom apeglm apeglm bbEstDisp
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_minimal labs coord_cartesian
#'
#' @export
estDisp <- function(sce, ct, pc = 2, genecluster) {
  if (missing(ct)) {
    cell_sum <- colSums(counts(sce)) %>% by(sce$x, FUN = mean)
    ct <- names(which.max(cell_sum))
  }

  if (missing(genecluster)) {
    cl <- metadata(sce)$geneCluster
    genecluster <- names(cl[which.max(cl)])
  }

  sce_sub <- sce[rowData(sce)$cluster == genecluster, sce$x == ct]
  x <- matrix(1, ncol = 1, nrow = dim(sce_sub)[2])
  theta.hat <- 100
  param <- cbind(theta.hat, counts(sce_sub))
  maxDisp <- 5000
  for (i in 1:5) {
    fit.mle <- apeglm(
      Y = assays(sce_sub)[["ase.mat"]], x = x, log.lik = NULL,
      param = param, no.shrink = TRUE, log.link = FALSE, method = "betabinC"
    )
    theta.hat <- bbEstDisp(
      success = assays(sce_sub)[["ase.mat"]],
      size = counts(sce_sub), x = x,
      beta = fit.mle$map, minDisp = .01, maxDisp = maxDisp
    )
  }
  gene_mean <- rowMeans(counts(sce_sub))
  est <- data.frame(mean = gene_mean, theta = theta.hat)
  est <- est[est$mean > 2 & est$theta < 100, ] # focus on genes with evidence of over-dispersion
  p <- ggplot(est, aes(mean, .data$theta)) +
    geom_point() +
    geom_smooth() +
    coord_cartesian(ylim = c(0, 100)) +
    theme_minimal() +
    labs(x = "gene mean", y = "theta")
  p
}
