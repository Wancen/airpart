#' @title Fit betabinomial on each cell type groups
#'
#' @description Fit betabinomial on each cell type groups
#'
#' @param sce A SingleCellExperiment containing assays
#' (\code{"ratio"}, \code{"counts"}) and colData (\code{"x"},
#' \code{"part"})
#' @param level the confidence interval required
#' @param ... Argument for the \code{\link[apeglm]{apeglm}}
#' functions.
#'
#' @return Posterior mean (\code{"ar"}) for allelic ratio estimator is returned in
#' rowData for each cell type as well as \code{"s"} value and confidence interval(\code{"lower"} and \code{"upper"}).
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = seq_len(4))
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' sce_sub <- allelicRatio(sce_sub)
#' @importFrom emdbook dbetabinom
#'
#' @export
allelicRatio <- function(sce, level = 0.95, ...) {
  nct <- nlevels(sce$x)
  x <- model.matrix(~x, colData(sce))
  theta.hat <- matrix(rep(100, dim(sce)[1])) # rough initial estimate of dispersion
  maxDisp <- 75
  niter <- 5
  for (i in seq_len(niter)) {
    param <- cbind(theta.hat[, 1], counts(sce))
    fit.mle <- apeglm(
      Y = assays(sce)[["a1"]], x = x, log.lik = NULL,
      param = param, no.shrink = TRUE, log.link = FALSE, method = "betabinC"
    )
    theta.hat <- bbEstDisp(
      success = assays(sce)[["a1"]],
      size = counts(sce), x = x,
      beta = fit.mle$map, minDisp = .01, maxDisp = maxDisp, se = TRUE
    )
  }
  theta <- theta.hat[, 1] # derive theta estimation given coefficient prior Cauchy(0,1) and betabinomial likelihood
  se <- theta.hat[, 2]
  ## approximate shrinkage from Efron and Morris
  log_disp <- log(theta)
  sigma_sampling2 <- mean(se^2)
  sigma_prior2 <- max((var(log_disp) - sigma_sampling2), 0)
  B <- as.numeric((sigma_sampling2) / (sigma_prior2 + sigma_sampling2))
  log_disp_0 <- mean(log_disp)
  theta2 <- exp((1 - B) * log_disp + B * log_disp_0)

  if ("coef" %in% colnames(colData(sce))) {
    ## use fused lasso estimates as prior
    coef <- unique(data.frame(x = sce$x, coef = sce$coef))$coef
  } else {
    ## use weighted mean as prior
    cl_ratio <- as.vector(unlist(assays(sce)[["ratio"]]))
    cl_total <- as.vector(unlist(assays(sce)[["counts"]]))
    dat <- data.frame(
      ratio = cl_ratio,
      part = factor(rep(sce$part, each = length(sce))),
      cts = cl_total
    )
    summary <- dat %>%
      group_by(.data$part) %>%
      summarise(coef = weighted.mean(.data$ratio, .data$cts, na.rm = TRUE)) %>%
      as.data.frame()
    coef <- merge(summary, metadata(sce)$partition, by = "part", sort = FALSE)$coef
    coef <- log(coef / (1 - coef))
  }

  param <- cbind(theta2, counts(sce))
  mean <- matrix(0, nrow = length(sce), ncol = nct) %>% `colnames<-`(paste("ar", levels(sce$x), sep = "_"))
  s <- matrix(0, nrow = length(sce), ncol = nct) %>% `colnames<-`(paste("svalue", levels(sce$x), sep = "_"))
  lower <- matrix(0, nrow = length(sce), ncol = nct) %>% `colnames<-`(paste("lower", paste0((1 - level) * 50, "pct"), levels(sce$x), sep = "_"))
  upper <- matrix(0, nrow = length(sce), ncol = nct) %>% `colnames<-`(paste("upper", paste0(100 - (1 - level) * 50, "pct"), levels(sce$x), sep = "_"))
  X <- model.matrix(~ x + 0, colData(sce))
  ## Shrink one cell type/time
  for (i in seq_len(nct)) {
    fit.mle <- apeglm(
      Y = assays(sce)[["a1"]], x = X, log.lik = betabinom.log.lik, coef = i,
      param = param, method = "general", interval.level = level,
      prior.control = list(
        no.shrink = setdiff(1:nct, i), prior.mean = coef, prior.scale = 1,
        prior.df = 1, prior.no.shrink.mean = 0, prior.no.shrink.scale = 15
      ), ...
    )
    mean[, i] <- (1 + exp(-fit.mle$map[, i]))^-1
    s[, i] <- fit.mle$svalue
    lower[, i] <- (1 + exp(-fit.mle$interval[, 1]))^-1
    upper[, i] <- (1 + exp(-fit.mle$interval[, 2]))^-1
  }
  est <- data.frame(round(mean, 3), format(s, digits = 3), round(lower, 3), round(upper, 3)) %>% `rownames<-`(rownames(sce))
  rowData(sce) <- merge(rowData(sce), est[, setdiff(colnames(est), colnames(rowData(sce)))], by = 0)[-1] %>%
    DataFrame() %>%
    `rownames<-`(rownames(sce))
  return(sce)
}

betabinom.log.lik <- function(y, x, beta, param, offset) {
  xbeta <- x %*% beta
  p.hat <- (1 + exp(-xbeta))^-1
  emdbook::dbetabinom(y, prob = p.hat, size = param[-1], theta = param[1], log = TRUE)
}
