#' Allelic ratio summary
#'
#' produce allelic ratio summary for each gene cluster
#'
#' @param sce SingleCellExperiment
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
#' @import magrittr
#' @importFrom dplyr summarise group_by
#'
#' @export
summaryAllelicRatio <- function(sce, genecluster) {
  if (missing(genecluster)) {
    genecluster <- seq_len(max(rowData(sce)$cluster))
  }
  res <- lapply(genecluster, function(i){
    sce_sub <- sce[rowData(sce)$cluster == i, ]
    cl_ratio <- as.vector(unlist(assays(sce_sub)[["ratio"]]))
    cl_total <- as.vector(unlist(assays(sce_sub)[["counts"]]))
    dat <- data.frame(ratio=cl_ratio,
                      x=factor(rep(sce_sub$x,each=length(sce_sub))),
                      cts=cl_total)
    summary <- dat %>%
      group_by(x) %>%
      summarise(weighted.mean=weighted.mean(ratio,cts,na.rm = TRUE),
                mean=mean(ratio,na.rm = TRUE),
                var=var(ratio,na.rm = TRUE)) %>%
      as.data.frame()
    summary
  })
  names(res) <- paste("gene cluster",genecluster,"with",metadata(sce)$geneCluster,"genes")
  metadata(sce)$summary <- res
  return(sce)
}
