#' Estimate overdispersion parameter of a betabinomial
#'
#' @param sce SingleCellExperiment with \code{a1} matrix and \code{counts}
#' @param pc pseudocount for calculating the smoothed
#' ratio in the preprocess step.
#' @param genecluster which gene cluster dispersion
#' parameter want to be estimated.
#' Default is the cluster with the most cells
#' @param type indicated whether to output
#' dispersion estimates as a plot or a value
#'
#' @return A ggplot object of the dispersion estimates over the mean
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = seq_len(4))
#' estDisp(sce)
#'
#' @importFrom apeglm apeglm bbEstDisp
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth
#' theme_minimal labs coord_cartesian
#' @importFrom stats model.matrix
#'
#' @export
estDisp <- function(sce, pc=2, genecluster, type = c("plot", "values")) {
  type <- match.arg(type, c("plot", "values"))
  if (missing(genecluster)) {
    cl <- metadata(sce)$geneCluster
    genecluster <- names(cl[which.max(cl)])
  }
  sce_sub <- sce[rowData(sce)$cluster == genecluster,]
  x <- model.matrix(~x, colData(sce))
  theta.hat <- 100
  param <- cbind(theta.hat, counts(sce_sub))
  maxDisp <- 5000
  for (i in seq_len(5)) {
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
  # focus on genes with evidence of over-dispersion
  est <- est[est$mean > 2 & est$theta < 100, ]

  # TODO: consider also outputting the dispersion estimates?
  # could do a `type="plot"` or `"values"` argument
  if (type == plot) {
    p <- ggplot(est, aes(mean, .data$theta)) +
      geom_point() +
      geom_smooth() +
      coord_cartesian(ylim = c(0, 100)) +
      theme_minimal() +
      labs(x = "gene mean", y = "theta")
    print(p)
  }

}
