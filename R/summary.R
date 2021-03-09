#' Allelic ratio summary
#'
#' produce allelic ratio summary for each gene cluster
#'
#' @param se SummarizedExpeirment
#' @param genecluster an optional vector of gene cluster IDs.
#' if nothing is given, all clusters' summary will be displayed
#'
#' @return a list of gene cluster's summary table contains
#' \itemize{
#'   \item{weighted.mean} {weighted mean of allelic ratio for the cell types}
#'   \item{mean} {mean allelic ratio for the cell types}
#'   \item{var} {variance of allelic ratio for the cell types}
#' }
#'
#' @examples summaryAllelicRatio(se,genecluster = c(1,3))
#'
#' @import tidyverse
#' @export
summaryAllelicRatio <- function(se,genecluster) {
  if (missing(genecluster)) {
    genecluster <- unique(rowData(se)$cluster)
  }
  res<-lapply(genecluster, function(i){
    se_sub<-se[rowData(se)$cluster == i, ]
    cl_ratio <- as.vector(unlist(assays(se_sub)[["ratio"]]))
    cl_total <- as.vector(unlist(assays(se_sub)[["total"]]))
    dat <- data.frame(ratio=cl_ratio,
                      x=factor(rep(colData(se_sub)$x,each=length(se_sub))),
                      cts=cl_total)
    summary <- dat %>%
      group_by(x) %>%
      summarise(weighted.mean=weighted.mean(ratio,cts,na.rm = T),
                mean=mean(ratio,na.rm = T),
                var=var(ratio,na.rm = T)) %>%
      as.data.frame()
    summary
  })
  names(res) <- paste("gene cluster",genecluster)
  metadata(se)$summary<-res
  return(se)
}



