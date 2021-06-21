#' Fit beta-binomial across cell types
#'
#' This function performs additional inference on the allelic ratio
#' across cell types, giving posterior mean and credible
#' intervals per cell type. A Cauchy prior is centered for
#' each cell type on the allelic ratio from the fused lasso
#' across all genes in the gene cluster
#' (or using a weighted means if the fused lasso is not provided).
#'
#' @param sce A SingleCellExperiment containing assays
#' (\code{"ratio"}, \code{"counts"}) and colData (\code{"x"},
#' \code{"part"})
#' @param formula The same \code{\link[stats]{formula}} object
#' used in \code{\link[airpart]{fusedLasso}}
#' @param level the level of credible interval (default is 0.95)
#' @param ... arguments to pass to \code{\link[apeglm]{apeglm}}
#' functions
#'
#' @return posterior mean (\code{"ar"}) for allelic ratio
#' estimate is returned in the rowData for each cell type,
#' as well as the \code{"s"} value and
#' credible interval (\code{"lower"} and \code{"upper"}).
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = seq_len(4))
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' sce_sub <- allelicRatio(sce_sub)
#' @importFrom emdbook dbetabinom
#' @importFrom stats as.formula
#'
#' @export
allelicRatio <- function(sce, formula, level = 0.95, ...) {
  if (missing(formula)) {
    formula <- ratio ~ p(x, pen = "gflasso")
  }
  add_covs <- grep("p\\(",
    attr(terms(formula), "term.labels"),
    invert = TRUE, value = TRUE
  )
  if (nlevels(sce$part) == 1) {
    x <- list(part = sce$part)
  } else {
    x <- list(part = model.matrix(~ part + 0, colData(sce)))
  }
  if (length(add_covs) > 0) {
    for (v in add_covs) {
      x[[v]] <- model.matrix(as.formula(paste("~", add_covs, "+0")), colData(sce))
      colnames(x[[v]]) <- levels(sce[[v]])
    }
  }
  x <- do.call(cbind, x)
  # rough initial estimate of dispersion
  theta.hat <- matrix(rep(100, dim(sce)[1]))
  maxDisp <- 75
  niter <- 5
  for (i in seq_len(niter)) {
    param <- cbind(theta.hat[, 1], counts(sce))
    fit.mle <- apeglm(
      Y = assays(sce)[["a1"]], x = x, log.lik = NULL,
      param = param, no.shrink = TRUE, log.link = FALSE,
      method = "betabinCR", interval.level = level,
    )
    theta.hat <- bbEstDisp(
      success = assays(sce)[["a1"]],
      size = counts(sce), x = x,
      beta = fit.mle$map, minDisp = .01, maxDisp = maxDisp, se = TRUE
    )
  }
  # construct offset if adjusted for batch, and sum them over different batch
  theta <- theta.hat[, 1]
  offset <- matrix(0, ncol = ncol(sce), nrow = nrow(sce))
  if (length(add_covs) > 0) {
    for (v in add_covs) {
      if (nrow(sce) == 1) {
        batch <- fit.mle$map[, match(levels(sce[[v]]), colnames(x))] %>%
          as.matrix() %>%
          t()
      } else {
        batch <- fit.mle$map[, match(levels(sce[[v]]), colnames(x))]
      }
      feat <- sce[[v]]
      offset <- offset + batch[, match(feat, colnames(batch))]
    }
  }
  if (nrow(sce) == 1) {
    param <- cbind(theta, counts(sce))
    coef <- 0
    res <- adp.shrink(sce, fit.mle, param, level, offset, coef, log.lik = NULL, method = "betabinCR", ...)
  } else {
    se <- theta.hat[, 2]
    ## approximate shrinkage using a formula in the style of
    # Efron and Morris's hierarchical model of normals
    log_disp <- log(theta)
    sigma_sampling2 <- mean(se^2)
    sigma_prior2 <- max((var(log_disp) - sigma_sampling2), 0)
    B <- as.numeric((sigma_sampling2) / (sigma_prior2 + sigma_sampling2))
    log_disp_0 <- mean(log_disp)
    theta2 <- exp((1 - B) * log_disp + B * log_disp_0)
    param <- cbind(theta2, counts(sce))
    ## estimating prior mean
    if ("coef" %in% colnames(colData(sce))) {
      ## use fused lasso estimates as prior
      coef <- unique(data.frame(x = sce$part, coef = sce$coef))$coef
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
        summarise(coef = weighted.mean(.data$ratio,
          .data$cts,
          na.rm = TRUE
        )) %>%
        as.data.frame()
      coef <- merge(summary, metadata(sce)$partition,
        by = "part", sort = FALSE
      )$coef
      coef <- log(coef / (1 - coef))
    }
    ## Shrink one cell type/time
    res <- adp.shrink(sce, fit.mle, param, level, offset, coef, log.lik = betabinom.log.lik, method = "general", ...)
  }
  mean <- res$mean
  s <- res$s
  lower <- res$lower
  upper <- res$upper
  est <- data.frame(
    round(mean, 3), format(s, digits = 3),
    round(lower, 3), round(upper, 3)
  ) %>%
    `rownames<-`(rownames(sce))
  est_sub <- est[, setdiff(colnames(est), colnames(rowData(sce)))]
  merged <- merge(rowData(sce), est_sub, by = 0) %>% DataFrame()
  order <- match(rownames(sce), merged$Row.names)
  rowData(sce) <- merged[order, ][-1] %>%
    `rownames<-`(rownames(sce))
  sce
}

## Betabinomial log link
betabinom.log.lik <- function(y, x, beta, param, offset) {
  xbeta <- x %*% beta + offset
  p.hat <- (1 + exp(-xbeta))^-1
  emdbook::dbetabinom(y,
    prob = p.hat, size = param[-1],
    theta = param[1], log = TRUE
  )
}

## derive adaptive shrinkage with prior
adp.shrink <- function(sce, fit.mle, param, level, offset, coef, log.lik, method, ...) {
  npart <- nlevels(sce$part)
  nct <- nlevels(sce$x)
  ## store results
  mean <- matrix(0, nrow = length(sce), ncol = nct) %>%
    `colnames<-`(paste("ar", levels(sce$x), sep = "_"))
  s <- matrix(0, nrow = length(sce), ncol = nct) %>%
    `colnames<-`(paste("svalue", levels(sce$x), sep = "_"))
  lower <- matrix(0, nrow = length(sce), ncol = nct) %>%
    `colnames<-`(paste("lower", paste0((1 - level) * 50, "pct"),
      levels(sce$x),
      sep = "_"
    ))
  upper <- matrix(0, nrow = length(sce), ncol = nct) %>%
    `colnames<-`(paste("upper", paste0(100 - (1 - level) * 50, "pct"),
      levels(sce$x),
      sep = "_"
    ))
  part <- unique(data.frame(x = sce$x, part = sce$part))$part
  if (nlevels(sce$part) == 1) {
    x <- data.frame(part = as.numeric(sce$part)) %>% as.matrix()
  } else {
    x <- model.matrix(~ part + 0, colData(sce))
  }
  for (i in seq_len(npart)) {
    mle <- cbind(fit.mle$map[, i], fit.mle$sd[, i])
    fit.post <- apeglm(
      Y = assays(sce)[["a1"]], x = x, param = param, log.lik = log.lik,
      coef = i, interval.level = level, offset = offset,
      log.link = FALSE, method = method,
      prior.control = list(
        no.shrink = setdiff(seq_len(ncol(x)), i), prior.mean = coef,
        prior.scale = 1, prior.df = 1,
        prior.no.shrink.mean = 0, prior.no.shrink.scale = 15
      ), ...
    )
    mean[, which(part == i)] <- (1 + exp(-fit.post$map[, i]))^-1
    s[, which(part == i)] <- fit.post$svalue
    lower[, which(part == i)] <- (1 + exp(-fit.post$interval[, 1]))^-1
    upper[, which(part == i)] <- (1 + exp(-fit.post$interval[, 2]))^-1
  }
  return(list(mean = mean, s = s, lower = lower, upper = upper))
}
