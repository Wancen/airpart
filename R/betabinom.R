#' @title Fit betabinomial on each cell type groups
#'
#' @description Fit betabinomial on each cell type groups
#'
#' @param se A SummarizedExpeirment containing assays (\code{"ratio"},
#' \code{"total"}) and colData (\code{"x"}, \code{"part"})
#' @param level the confidence interval required
#'
#' @importFrom  VGAM vglm betabinomial Coef confintvglm
#'
#' @export
betaBinom<-function(se,level=0.95,...){
  cl_ratio <- as.vector(unlist(assays(se)[["ratio"]]))
  cl_total <- as.vector(unlist(assays(se)[["total"]]))
  dat <- data.frame(ratio=cl_ratio,
                    x=factor(rep(se$x,each=length(se))),
                    cts=cl_total,
                    part=rep(se$part,each=length(se)))
  # Modeling each group separately because they may have different scale of over-dispersion
  estimator<-sapply (1:nlevels(se$part), function(m){
    suppressWarnings(bb<-VGAM::vglm(cbind(ratio*cts, cts-ratio*cts) ~1, VGAM::betabinomial, data = dat[which(dat$part==m),],
                   trace = F))
    coef_bb<-VGAM::Coef(bb)[-2] # betabinomial estimator
    rho<-VGAM::Coef(bb)[2]
    suppressWarnings(confint_bb<-VGAM::confintvglm(bb,matrix=T,level = level)[-2,]) #ci
    confint_wilcoxon<-1/(1+exp(-confint_bb))
    return(list(coef_bb,confint_wilcoxon))
  })
  coef <- as.vector(do.call(rbind, estimator[seq(1,length(estimator), by = 2)]))
  est0 <- data.frame(part=factor(seq_len(length(coef))), estimator=round(coef,3)) # allelic ratio estimator

  confint <- matrix(do.call(rbind, estimator[seq(2,length(estimator), by = 2)]),ncol=2) %>%
    as.data.frame() %>%
    setNames(names(estimator[[2]]))
  ci <- cbind(part=factor(seq_len(length(coef))),round(confint,3))
  coldata<-Reduce(function(x,y) merge(x = x, y = y, by = "part"),
         list(metadata(se)$partition, est0, ci))
  # colData(se) <-coldata %>% DataFrame() %>% setNames(colnames(coldata)) # combine with partition label
  metadata(se)$estimator<-coldata
  return(se)
}
