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
  cl <- data.frame(metadata(sce)$partition) %>% select(-.data$x)
  consens_part <- cl_consensus(
    cl_ensemble(list=apply(cl,2,as.cl_hard_partition)),
    method = "SM"
  )
  class <- max.col(consens_part$.Data)
  class <- match(class,unique(class))
  cl <- data.frame(part=factor(class), x=levels(sce$x))
  colData(sce) <- merge(colData(sce),cl,by="x")%>% DataFrame() %>% `row.names<-`(colnames(sce))
  metadata(sce)$partition <- cl
  return(sce)
}
