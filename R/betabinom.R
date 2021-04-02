#' @title Fit betabinomial on each cell type groups
#'
#' @description Fit betabinomial on each cell type groups
#'
#' @param sce A SingleCellExperiment containing assays (\code{"ratio"},
#' \code{"counts"}) and colData (\code{"x"}, \code{"part"})
#' @param level the confidence interval required
#' @param method the method to used for construct confidence interval. The default is \code{"normal"}.
#' \code{"bootstrap"} is likely to be more accurate especially when the data is highly overdispersed,
#' however it is computationally intensive.
#' @param type the type of intervals required for bootstrap. The values should be one of the values
#' \code{c("norm","basic","perc")}
#' @param R The number of bootstrap replicates.
#' @param trace logical indicating if output should be produced for each iteration.
#' @param ... Argument for the \code{\link[boot]{boot}} functions.
#'
#' @return A matrix allelic ratio estimator is returned in
#' metadata \code{"estimator"} which includes \code{"estimate"}, \code{"confidence interval"}
#' and \code{"std.error"} if use normal approximation.
#'
#' @examples
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = seq_len(4))
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#'
#' # Use normal approximation to calculate Wald-type confidence intervals
#' sce_sub <- allelicRatio(sce_sub)
#'
#' # Alternative with bootstrap
#' sce_sub <- allelicRatio(sce_sub,
#'   method = "bootstrap", R = 5, type="norm",
#'   parallel = "multicore", ncpus = 4
#' )
#'
#' @importFrom boot boot boot.ci
#' @importFrom stats setNames coef pt qnorm
#' @importFrom VGAM vglm betabinomial confintvglm
#' @importFrom matrixStats colSds
#'
#' @export
allelicRatio <- function(sce, level = 0.95, method = c("normal", "bootstrap"), type = "perc",
                         R, trace = TRUE,...) {
  method <- match.arg(method, c("normal", "bootstrap"))[1]
  cl_ratio <- as.vector(unlist(assays(sce)[["ratio"]]))
  cl_total <- as.vector(unlist(counts(sce)))

  stopifnot(c("x","part") %in% names(colData(sce)))

  dat <- data.frame(
    ratio = cl_ratio,
    x = factor(rep(sce$x, each = length(sce))),
    cts = cl_total,
    part = rep(sce$part, each = length(sce))
  )
  dat <- dat[!is.nan(dat$ratio), ]
  n <- table(dat$x)
  if (method == "bootstrap") {
    if (missing(R)) {
      stop("No number of bootstrap replicates")
    }
    boot <- boot(dat, statistic = boot_ci, R = R, strata = dat$part, trace = trace, ...)
    confint <- vapply(seq_len(nlevels(dat$x)), function(m) {
      boot_ci <- boot.ci(boot, type = type, index = m, conf = level)
      ci <- boot_ci[[length(boot_ci)]]
      return(ci[(length(ci) - 1):length(ci)])
    }, double(2))
    statistics <- (as.vector(boot$t0) - 0.5) / colSds(boot$t)
    pvalue <- data.frame(pvalue = format(2 * pt(abs(statistics), df = n - 1, lower.tail = FALSE), digits = 4))
    coef <- cbind(estimate = as.vector(boot$t0), std.error = colSds(boot$t)) %>% as.data.frame()
    confint <- t(confint) %>%
      as.data.frame() %>%
      setNames(paste(c((1 - level) * 50, 100 - (1 - level) * 50), "%"))
  } else if (method == "normal") {
    estimator <- betaBinom(dat, ci = TRUE, level = level, trace = trace)
    coef <- do.call(c, lapply(estimator, `[[`, 1))
    confint <- matrix(do.call(rbind, lapply(estimator, `[[`, 2)), ncol = 2) %>%
      as.data.frame() %>%
      setNames(paste(c((1 - level) * 50, 100 - (1 - level) * 50), "%"))
    aa <- (1 - level) / 2
    fac <- qnorm(aa)
    se <- (coef - confint[, 1]) / abs(fac)
    statistics <- (coef - 0.5) / se
    coef <- cbind(estimate = coef, std.error = se) %>% as.data.frame()
    pvalue <- data.frame(pvalue = format(2 * pt(abs(statistics), df = n - 1, lower.tail = FALSE), digits = 4))
  }
  order <- dat %>%
    group_by(.data$part) %>%
    summarise(x = unique(.data$x))
  est <- cbind(x = order$x, round(coef, 3)) # allelic ratio estimator
  pvalue <- cbind(x = order$x, pvalue = pvalue)
  ci <- cbind(x = order$x, round(confint, 3))
  coldata <- Reduce(
    function(x, y) merge(x = x, y = y, by = "x"),
    list(metadata(sce)$partition, est, pvalue, ci)
  )
  coldata <- coldata[order(match(coldata$x, levels(sce$x))), ] # change cell type order
  metadata(sce)$estimator <- coldata
  return(sce)
}

# betabinomial bootstrap estimator
betaBinom <- function(data, ci, level, trace) {
  # modeling each group separately because they may
  # have different scale of over-dispersion
  res <- lapply(seq_len(nlevels(data$part)), function(m) {
    data2 <- data[which(data$part == m), ]
    if (length(unique(data2$x)) == 1) {
      suppressWarnings({
        fit <- VGAM::vglm(cbind(ratio * cts, cts - ratio * cts) ~ 1,
                 VGAM::betabinomial(lmu = "identitylink", lrho = "identitylink"),
                 data = data2, trace = trace)
      })
    } else {
      suppressWarnings({
        fit <- VGAM::vglm(cbind(ratio * cts, cts - ratio * cts) ~ x,
                 VGAM::betabinomial(lmu = "identitylink", lrho = "identitylink"),
                 data = data2, trace = trace)
      })
    }
    coef <- coef(fit)[-2]
    coef <- coef + c(0, rep(coef[1], (length(coef) - 1)))
    if (ci) {
      suppressWarnings(confint_bb <- VGAM::confintvglm(fit, matrix = TRUE, level = level)[-2, ]) # ci
      confint <- confint_bb + matrix(c(0, 0, rep(coef[1], 2 * (length(coef) - 1))), byrow = TRUE, ncol = 2)
      return(list(coef, confint))
    } else {
      return(unname(coef))
    }
  })

  # return a matrix of coefficients (columns) if the CI were not requested
  if (!ci) {
    res <- do.call(c, res)
  }

  return(res)
}

# boostrap helper function
boot_ci <- function(data, indices, trace) {
  data_b <- data[indices, ]
  coef <- betaBinom(data_b, ci = FALSE, trace=trace)
  return(unname(coef))
}
