#' @import  VGAM
#' @export
#' @title Fit betabinomial on each cell type groups
#'
#' @description Fit betabinomial on each cell type groups
#'
#' @param data A data frame containing:
#' the allelic ratio (\code{"ratio"}),
#' and the total counts for weighting (\code{"cts"})
#' the cell type group (\code{"part"})
#' @param part the column name of cell type group
#'
#' @return a list with the following elements:
#' \item{coef_bb}{the allelic ratio estimator for each cell type group}
#' \item{confint_wilcoxon}{the allelic ratio confindence interval for each cell type group}
#' \item{rho}{the correlation between the N individuals within a group. rho is given by 1/(1 + alpha + beta) and is
#' It is known as the over-dispersion parameter}
betabinom<-function(data,part,...){
  cl<-data[[part]]
  # Modeling each group separately because they may have different scale of over-dispersion
  estimator<-sapply (1:max(cl), function(m){
    bb<-VGAM::vglm(cbind(ratio*cts, cts-ratio*cts) ~1, VGAM::betabinomial, data = data[which(cl==m),], trace = F)
    coef_bb<-VGAM::Coef(bb)[-2] # betabinomial estimator
    rho<-VGAM::Coef(bb)[2]
    confint_bb<-VGAM::confintvglm(bb,matrix=T)[-2,] #ci
    confint_wilcoxon<-1/(1+exp(-confint_bb))
    return(list(coef_bb,confint_wilcoxon,rho))
  })
}
