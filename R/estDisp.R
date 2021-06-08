#' Estimate overdispersion parameter of a beta-binomial
#'
#' @param sce SingleCellExperiment with \code{a1} matrix and \code{counts}
#' @param genecluster the gene cluster for which to estimate the
#' over-dispersion parameter. Default is the cluster with the most cells
#' @param type whether to output the over-dispersion estimates as a plot or a value
#'
#' @return A ggplot object of the dispersion estimates over the mean,
#' or a data.frame of the mean and dispersion estimates (theta)
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = seq_len(4))
#' estDisp(sce, genecluster = 1)
#' @importFrom apeglm apeglm bbEstDisp
#' @importFrom stats model.matrix
#'
#' @export
estDisp <- function(sce, genecluster, type = c("plot", "values")) {
  type <- match.arg(type, c("plot", "values"))
  stopifnot("cluster" %in% names(rowData(sce)))
  if (missing(genecluster)) {
    cl <- metadata(sce)$geneCluster
    genecluster <- names(cl[which.max(cl)])
  }
  sce_sub <- sce[rowData(sce)$cluster == genecluster, ]
  x <- model.matrix(~x, colData(sce))
  theta.hat <- 100
  param <- cbind(theta.hat, counts(sce_sub))
  maxDisp <- 5000
  for (i in seq_len(5)) {
    param <- cbind(theta.hat, counts(sce_sub))
    fit.mle <- apeglm(
      Y = assays(sce_sub)[["a1"]], x = x, log.lik = NULL,
      param = param, no.shrink = TRUE, log.link = FALSE, method = "betabinC"
    )
    theta.hat <- bbEstDisp(
      success = assays(sce_sub)[["a1"]],
      size = counts(sce_sub), x = x,
      beta = fit.mle$map, minDisp = .01, maxDisp = maxDisp
    )
  }
  gene_mean <- rowMeans(counts(sce_sub))
  est <- data.frame(mean = gene_mean, theta = theta.hat)
  ## focus on genes with evidence of over-dispersion
  est <- est[est$mean > 2 & est$theta < 100, ]
  if (type == "plot") {
    p <- ggplot(est, aes(mean, .data$theta)) +
      geom_point() +
      geom_smooth() +
      coord_cartesian(ylim = c(0, 100)) +
      theme_minimal() +
      labs(x = "gene mean", y = "theta")
    p
  } else {
    est
  }
}
