#' @title Fit betabinomial on each cell type groups
#'
#' @description Fit betabinomial on each cell type groups
#'
#' @param sce A SingleCellExperiment containing assays (\code{"ratio"},
#' \code{"counts"}) and colData (\code{"x"}, \code{"part"})
#' @param level the confidence interval required
#' @param method the method to used for construct confidence interval. The default is the first method.
#' "bootstrap" is likely to be more accurate especially when the data is highly overdispersed,
#' however it is computationally intensive.
#' @param R The number of bootstrap replicates.
#' @param ... Argument for the \code{\link[boot]{boot}} functions.
#'
#' @return A matrix allelic ratio estimator is returned in
#' metadata \code{"estimator"}
#'
#' @examples
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#'
#' # Use normal approximation to calculate Wald-type confidence intervals
#' sce_sub <- allelicRatio(sce_sub)
#'
#' # Alternative with bootstrap
#' \dontrun{
#' sce_sub <- allelicRatio(sce_sub,
#'   method = "bootstrap", R = 200,
#'   parallel = "multicore", ncpus = 4
#' )
#' }
#'
#' @importFrom VGAM vglm betabinomial Coef confintvglm
#' @importFrom boot boot boot.ci
#' @importFrom stats setNames
#'
#' @export
allelicRatio <- function(sce, level = 0.95, method = c("Normal", "bootstrap"),
                         R, ...) {
  method <- match.arg(method, c("Normal", "bootstrap"))[1]
  cl_ratio <- as.vector(unlist(assays(sce)[["ratio"]]))
  cl_total <- as.vector(unlist(counts(sce)))
  dat <- data.frame(
    ratio = cl_ratio,
    x = factor(rep(sce$x, each = length(sce))),
    cts = cl_total,
    part = rep(sce$part, each = length(sce))
  )
  if (method == "bootstrap") {
    boot <- boot(dat, statistic = boot_ci, R = R, strata = dat$part, level = level, ...)
    confint <- sapply(1:nlevels(dat$part), function(m) {
      boot_ci <- boot.ci(boot, type = "perc", index = m, conf = level)
      ci <- boot_ci[[length(boot_ci)]]
      return(ci[(length(ci) - 1):length(ci)])
    })
    coef <- boot$t0
    confint <- t(confint) %>%
      as.data.frame() %>%
      setNames(paste(c((1 - level) * 50, 100 - (1 - level) * 50), "%"))
  } else {
    estimator <- betaBinom(dat, level = level)
    coef <- as.vector(do.call(rbind, estimator[seq(1, length(estimator), by = 2)]))
    confint <- matrix(do.call(rbind, estimator[seq(2, length(estimator), by = 2)]), ncol = 2) %>%
      as.data.frame() %>%
      setNames(names(estimator[[2]]))
  }
  est <- data.frame(part = factor(seq_len(length(coef))), estimator = round(coef, 3)) # allelic ratio estimator
  ci <- cbind(part = factor(seq_len(length(coef))), round(confint, 3))
  coldata <- Reduce(
    function(x, y) merge(x = x, y = y, by = "part"),
    list(metadata(sce)$partition, est, ci)
  )
  # colData(sce) <-coldata %>% DataFrame() %>% setNames(colnames(coldata)) # combine with partition label
  metadata(sce)$estimator <- coldata
  return(sce)
}

betaBinom <- function(data, ci = TRUE, level) {
  # Modeling each group separately because they may have different scale of over-dispersion
  res <- sapply(1:nlevels(data$part), function(m) {
    suppressWarnings(bb <- VGAM::vglm(cbind(ratio * cts, cts - ratio * cts) ~ 1, VGAM::betabinomial,
      data = data[which(data$part == m), ],
      trace = F
    ))
    coef_bb <- VGAM::Coef(bb)[-2] # betabinomial estimator
    if (ci) {
      suppressWarnings(confint_bb <- VGAM::confintvglm(bb, matrix = T, level = level)[-2, ]) # ci
      confint_wilcoxon <- 1 / (1 + exp(-confint_bb))
      return(list(coef_bb, confint_wilcoxon))
    } else {
      return(unname(coef_bb))
    }
  })
  return(res)
}

boot_ci <- function(data, indices, level = level) {
  data_b <- data[indices, ]
  coef <- betaBinom(data_b, ci = FALSE, level)
  return(coef)
}
