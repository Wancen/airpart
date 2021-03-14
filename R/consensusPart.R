#' Consensus Partitions
#'
#' Derive consensus partitions of an ensemble fused lasso partitions.
#'
#' @param sce SingleCellExperiment
#'
#' @importFrom clue cl_consensus cl_ensemble as.cl_hard_partition
#' @importFrom dplyr select
#'
#' @export
consensusPart <- function(sce) {
  cl <- data.frame(metadata(sce)$partition) %>% select(-x)
  consens_part <- cl_consensus(
    cl_ensemble(list=apply(cl,2,as.cl_hard_partition)),
    method = "SM"
  )
  class <- max.col(consens_part$.Data)
  class <- match(class,unique(class))
  cl <- data.frame(part=factor(class), x=levels(sce_sub$x))
  colData(sce) <- merge(colData(sce),cl,by="x")%>% DataFrame() %>% `row.names<-`(colnames(sce_sub))
  metadata(sce)$partition <- cl
  return(sce)
}
