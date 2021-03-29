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
#' @importFrom gamlss gamlss
#' @importFrom gamlss.dist BB
#' @importFrom broom.mixed tidy
#' @importFrom boot boot boot.ci
#' @importFrom stats setNames coef
#' @importFrom VGAM vglm betabinomial
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
  dat <- dat[!is.nan(dat$ratio), ]
  if (method == "bootstrap") {
    boot <- boot(dat, statistic = boot_ci, R = R, strata = dat$part, ...)
    confint <- sapply(1:nlevels(dat$x), function(m) {
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
    coef <- estimator[,c("estimate","std.error")]
    confint <- estimator[,c("conf.low","conf.high")] %>%
      setNames(paste(c((1 - level) * 50, 100 - (1 - level) * 50), "%"))
  }
  est <- cbind(x = metadata(sce)$partition$x, round(coef, 3)) # allelic ratio estimator
  ci <- cbind(x = metadata(sce)$partition$x, round(confint, 3))
  coldata <- Reduce(
    function(x, y) merge(x = x, y = y, by = "x"),
    list(metadata(sce)$partition, est, ci)
  )
  coldata <- coldata[order(match(coldata$x, levels(sce$x))),] # change cell type order
  metadata(sce)$estimator <- coldata
  return(sce)
}

# betabinomial estimator
betaBinom <- function(data,level) {
  # Modeling each group separately because they may have different scale of over-dispersion
  suppressWarnings(fit <- gamlss::gamlss(cbind(ratio * cts, cts - ratio * cts) ~ x+0,
                                         sigma.formula = ~part+0, # allow different dispersion for each group
                                         data = data,
                                         family = gamlss.dist::BB(mu.link = "identity")))
  td <- broom.mixed::tidy(fit,conf.int = T,conf.level = level)
  return(td[which(td$parameter=="mu"),c(".rownames","estimate","std.error","conf.low","conf.high")])
}

# betabinomial bootstrap estimator
betaBinom_boot <- function(data) {
  # Modeling each group separately because they may have different scale of over-dispersion
  res <- sapply(1:nlevels(data$part), function(m) {
    data2<-data[which(data$part == m), ]
    if(length(unique(data2$x))==1){
      suppressWarnings(fit <- VGAM::vglm(cbind(ratio * cts, cts - ratio * cts) ~ 1, VGAM::betabinomial,
                                         data = data2,
                                         trace = F
      ))
    }else{
      suppressWarnings(fit <- VGAM::vglm(cbind(ratio * cts, cts - ratio * cts) ~ x, VGAM::betabinomial,
                                         data = data2,
                                         trace = F
      ))
    }
    coef<-coef(fit)[-2]
    coef<-coef+c(0,rep(coef[1],(length(coef)-1)))
    coef <- 1 / (1 + exp(-coef))# betabinomial estimator
    names(coef) <- unique(data2$x)
    return(coef)
  })
  coef <- unlist(res)
  coef <- coef[order(match(names(coef), levels(data$x)))]
  return(unname(coef))
}
# boostrap helper function
boot_ci <- function(data, indices) {
  data_b <- data[indices, ]
  coef <- betaBinom_boot(data_b)
  return(unname(unlist(coef)))
}
