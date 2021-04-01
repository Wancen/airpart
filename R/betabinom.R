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
#' @importFrom broom.mixed tidy
#' @importFrom boot boot boot.ci
#' @importFrom stats setNames coef pt qnorm
#' @importFrom VGAM vglm betabinomial confintvglm
#' @importFrom matrixStats colSds
#'
#' @export
allelicRatio <- function(sce, level = 0.95, method = c("Normal", "bootstrap"),
                         R, trace = TRUE, ...) {
  method <- match.arg(method, c("Normal", "bootstrap"))[1]
  cl_ratio <- as.vector(unlist(assays(sce)[["ratio"]]))
  cl_total <- as.vector(unlist(counts(sce)))
  dat <- data.frame(
    ratio = cl_ratio,
    x = factor(rep(sce$x, each = length(sce))),
    cts = cl_total,
    part = rep(sce$part, each = length(sce))
  )
  dat <- dat[!is.nan(dat$ratio), ]
  n <- table(dat$x)
  if (method == "bootstrap") {
    boot <- boot(dat, statistic = boot_ci, R = R, strata = dat$part,trace=trace, ...)
    confint <- sapply(1:nlevels(dat$x), function(m) {
      boot_ci <- boot.ci(boot, type = "perc", index = m, conf = level)
      ci <- boot_ci[[length(boot_ci)]]
      return(ci[(length(ci) - 1):length(ci)])
    })
    statistics <- (boot$t0-0.5)/colSds(boot$t)
    pvalue <- data.frame(pvalue = format(2 * pt(abs(statistics), df=n-1,lower.tail = FALSE),digits = 4))
    coef <- cbind(estimate=boot$t0,std.error=colSds(boot$t)) %>% as.data.frame()
    confint <- t(confint) %>%
      as.data.frame() %>%
      setNames(paste(c((1 - level) * 50, 100 - (1 - level) * 50), "%"))
  } else {
    estimator <- betaBinom(dat, ci = TRUE, level = level, trace=trace)
    coef <- unlist(estimator[seq(1, length(estimator), by = 2)])
    confint <- matrix(do.call(rbind, estimator[seq(2, length(estimator), by = 2)]), ncol = 2) %>%
      as.data.frame() %>%
      setNames(paste(c((1 - level) * 50, 100 - (1 - level) * 50), "%"))
    aa <- (1 - level) / 2
    fac <- qnorm(aa)
    se <- (coef - confint[,1])/abs(fac)
    statistics <- (coef-0.5)/se
    coef <- cbind(estimate=coef,std.error=se) %>% as.data.frame()
    pvalue <- data.frame(pvalue = format(2 * pt(abs(statistics), df=n-1,lower.tail = FALSE),digits = 4))
  }
  order <- dat %>% group_by(.data$part) %>% summarise(x=unique(.data$x))
  est <- cbind(x = order$x, round(coef, 3)) # allelic ratio estimator
  pvalue <- cbind(x = order$x, pvalue= pvalue)
  ci <- cbind(x = order$x, round(confint, 3))
  coldata <- Reduce(
    function(x, y) merge(x = x, y = y, by = "x"),
    list(metadata(sce)$partition, est, pvalue, ci)
  )
  coldata <- coldata[order(match(coldata$x, levels(sce$x))),] # change cell type order
  metadata(sce)$estimator <- coldata
  return(sce)
}

# betabinomial estimator
# betaBinom <- function(data,ci,level) {
#   # Modeling each group separately because they may have different scale of over-dispersion
  # suppressWarnings(fit <- gamlss::gamlss(cbind(ratio * cts, cts - ratio * cts) ~ x+0,
  #                                        sigma.formula = ~part+0, # allow different dispersion for each group
  #                                        data = data,
  #                                        family = gamlss.dist::BB("identity")))
# #   try <- tryCatch(
# #     {
#       td <- broom.mixed::tidy(fit,conf.int = T,conf.level = level)
#       TRUE
#     },
#     error = function(e) {
#         p <- fit$mu.df
#         df_r = fit$df.residual
#         w = fit$weights
#         idx = w > 0
#         dispersion = sum(w[idx] * fit$residuals[idx]^2)/df_r
#         p1 = seq_len(p)
#         nm <- names(fit$mu.coefficients[fit$mu.qr$pivot[p1]])
#         covmat = dispersion * chol2inv(fit$mu.qr$qr[p1, p1])
#         dimnames(covmat) = list(nm, nm)
#         se <- sqrt(diag(covmat))
#         cf <- fit$mu.coefficients
#         aa <- (1 - level) / 2
#         aa <- c(aa, 1 - aa)
#         fac <- qnorm(aa)
#         ci <- cf + se %o% fac
#         td <- data.frame(`.rownames`=nm,estimate=cf,`std.error`=se,`conf.low`= ci[,1],`conf.high`=ci[,2])
#     }
#   )
#   td$estimate <- 1 / (1 + exp(-td$estimate))
#   td$conf.low <- 1 / (1 + exp(-td$conf.low))
#   td$conf.high <- 1 / (1 + exp(-td$conf.high))
#   if(ci){
#     return(td[which(td$parameter=="mu"),c(".rownames","estimate","std.error","conf.low","conf.high")])
#   }else{
#     return(unname(unlist(td[which(td$parameter=="mu"),c("estimate")])))
#   }
# }

# betabinomial bootstrap estimator
betaBinom <- function(data, ci, level,trace) {
  # Modeling each group separately because they may have different scale of over-dispersion
  res <- sapply(1:nlevels(data$part), function(m) {
    data2<-data[which(data$part == m), ]
    if(length(unique(data2$x))==1){
      suppressWarnings(fit <- VGAM::vglm(cbind(ratio * cts, cts - ratio * cts) ~ 1, VGAM::betabinomial(lmu = "identitylink",
                                                                                                       lrho = "identitylink"),
                                         data = data2,
                                         trace = trace
      ))
    }else{
      suppressWarnings(fit <- VGAM::vglm(cbind(ratio * cts, cts - ratio * cts) ~ x, VGAM::betabinomial(lmu = "identitylink",
                                                                                                       lrho = "identitylink"),
                                         data = data2,
                                         trace = trace
      ))
    }
    coef<-coef(fit)[-2]
    coef<-coef+c(0,rep(coef[1],(length(coef)-1)))
    if (ci) {
      suppressWarnings(confint_bb <- VGAM::confintvglm(fit, matrix = T, level = level)[-2, ]) # ci
      confint <- confint_bb + matrix(c(0,0,rep(coef[1],2*(length(coef)-1))),byrow = T,ncol = 2)
      return(list(coef, confint))
    } else {
      return(unname(coef))
    }
  })
  return(res)
}
# boostrap helper function
boot_ci <- function(data, indices,trace) {
  data_b <- data[indices, ]
  coef <- betaBinom(data_b,ci=FALSE,trace)
  return(unname(unlist(coef)))
}
