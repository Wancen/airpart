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
#' @param nogroup Indicate whether there is previous group step. Default is FALSE represents either
#' Generalized fused lasso or nonparametric method is used to derive partition. nogroup == TRUE means
#' partition step is not needed to estimate allelic ratio.
#' @param level the level of credible interval (default is 0.95)
#' @param DAItest Indicate whether to do Likelihood Ratio test on differential allelic imbalance(DAI)
#' or equivalent to heterogeneity.
#' @param ... arguments to pass to \code{\link[apeglm]{apeglm}}
#' functions
#'
#' @return posterior mean (\code{"ar"}) for allelic ratio
#' estimate is returned in the rowData for each cell type,
#' as well as the \code{"s"} value, \code{"fsr"} false sign rate and
#' credible interval (\code{"lower"} and \code{"upper"}). One can use \code{"fsr"} < 0.005 or
#' credible intervals contain 0.5 or not for AI test significance. \code{"p.value"} shows DAI test result
#' and \code{"adj.p.value"} is false discovery rate corrected p values.
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = seq_len(4))
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' sce_sub <- allelicRatio(sce_sub, DAItest = TRUE)
#' @importFrom emdbook dbetabinom
#' @importFrom stats as.formula
#'
#' @export
allelicRatio <- function(sce, formula, nogroup = FALSE, level = 0.95, DAItest = FALSE, ...) {
  if (missing(formula)) {
    formula <- ratio ~ p(x, pen = "gflasso")
  }
  add_covs <- grep("p\\(",
                   attr(terms(formula), "term.labels"),
                   invert = TRUE, value = TRUE
  )
  if (isTRUE(nogroup)){
    sce$part = sce$x
  }

  # derive alternative hypothesis result
  x <- designMatrix(sce, add_covs)
  res <- part(sce, x, add_covs, level = level)

  if(isTRUE(DAItest)){
    npart = nlevels(sce$part)
    # reconstruct coefficient matrix
    coef = res$beta[,!duplicated(res$beta[1,])]
    ratio <- inv.logit(coef %*% t(x) + res$offset)
    loglik1 <- rowSums(emdbook::dbetabinom(assays(sce)[["a1"]], ratio, assays(sce)[["counts"]], res$theta, log=TRUE))
    # derive null hypothesis result
    sce0 <- sce
    sce0$part = factor(1)
    x0 <- designMatrix(sce0, add_covs)
    res0 <- part(sce0, x0, add_covs, level = level)
    coef0 = res0$beta[,!duplicated(res0$beta[1,])]
    ratio0 <- inv.logit(coef0 %*% t(x0) + res0$offset)
    loglik0 <- rowSums(emdbook::dbetabinom(assays(sce0)[["a1"]], ratio0, assays(sce0)[["counts"]], res0$theta, log=TRUE))
    # Construct LR
    LR = -2 * (loglik0-loglik1)
    rowData(sce)$p.value = pchisq(LR, npart-1, lower.tail = FALSE)
    rowData(sce)$adj.p.value = p.adjust(rowData(sce)$p.value, method = "fdr")
  }
  mean <- inv.logit(res$beta)
  s <- res$s
  fsr <- res$fsr
  lower <- res$lower
  upper <- res$upper
  est <- data.frame(
    round(mean, 3), format(s, digits = 3), format(fsr, digits = 3),
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

## construct design matrix
designMatrix <- function(sce, add_covs){
  if (nlevels(sce$part) == 1) {
    X <- list(part = sce$part)
  } else {
    X <- list(part = model.matrix(~ part + 0, colData(sce)))
  }

  if (length(add_covs) > 0) {
    for (v in add_covs) {
      X[[v]] <- model.matrix(as.formula(paste("~", add_covs, "+0")), colData(sce))
      colnames(X[[v]]) <- levels(sce[[v]])
    }
  }
  X <- do.call(cbind, X)
  return(X)
}


part <-function(sce, x, add_covs, level, ...){
  # rough initial estimate of dispersion
  theta.hat <- matrix(rep(100, dim(sce)[1]))
  maxDisp <- 75
  niter <- 5
  for (i in seq_len(niter)) {
    param <- cbind(theta.hat[, 1], assays(sce)[["counts"]])
    fit.mle <- apeglm(
      Y = assays(sce)[["a1"]], x = x, log.lik = NULL,
      param = param, no.shrink = TRUE, log.link = FALSE,
      method = "betabinCR", interval.level = level,
    )
    theta.hat <- bbEstDisp(
      success = assays(sce)[["a1"]],
      size = assays(sce)[["counts"]], x = x,
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
    param <- cbind(theta, assays(sce)[["counts"]])
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
    param <- cbind(theta2, assays(sce)[["counts"]])
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
      coef <- merge(summary, unique(data.frame(part = sce$part, x = sce$x)),
                    by = "part", sort = FALSE
      )$coef
      coef <- log(coef / (1 - coef))
    }
    ## Shrink one cell type/time
    res <- adp.shrink(sce, fit.mle, param, level, offset, coef, log.lik = betabinom.log.lik, method = "general", ...)
  }
  return(res)
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
  beta <- matrix(0, nrow = length(sce), ncol = nct) %>%
    `colnames<-`(paste("ar", levels(sce$x), sep = "_"))
  s <- matrix(0, nrow = length(sce), ncol = nct) %>%
    `colnames<-`(paste("svalue", levels(sce$x), sep = "_"))
  fsr <- matrix(0, nrow = length(sce), ncol = nct) %>%
    `colnames<-`(paste("fsr", levels(sce$x), sep = "_"))
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
  for (j in seq_len(npart)) {
    mle <- cbind(fit.mle$map[, j], fit.mle$sd[, j])
    fit.post <- apeglm(
      Y = assays(sce)[["a1"]], x = x, param = param, log.lik = log.lik,
      coef = j, interval.level = level, offset = offset,
      log.link = FALSE, method = method,
      prior.control = list(
        no.shrink = setdiff(seq_len(ncol(x)), j), prior.mean = coef,
        prior.scale = 1, prior.df = 1,
        prior.no.shrink.mean = 0, prior.no.shrink.scale = 15
      ), ...
    )
    beta[, which(part == levels(part)[j])] <- fit.post$map[, j]
    s[, which(part == levels(part)[j])] <- fit.post$svalue
    fsr[, which(part == levels(part)[j])] <- fit.post$fsr
    lower[, which(part == levels(part)[j])] <- inv.logit(fit.post$interval[, 1])
    upper[, which(part == levels(part)[j])] <- inv.logit(fit.post$interval[, 2])
  }
  return(list(beta = beta, s = s, fsr = fsr, lower = lower, upper = upper, offset = offset, theta = param[,1]))
}

inv.logit <- function(x) (1 + exp(-x))^-1
