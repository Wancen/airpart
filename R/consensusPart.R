#' Consensus Partitions
#'
#' Derive consensus partitions of an ensemble fused lasso partitions.
#'
#' @param se SummarizedExpeirment
#'
#' @importFrom clue cl_consensus cl_ensemble as.cl_hard_partition
#' @importFrom dplyr select
#' 
#' @export
consensusPart <- function(se) {
  cl <- data.frame(metadata(se)$partition) %>% select(-x)
  consens_part <- cl_consensus(
    cl_ensemble(list=apply(cl,2,as.cl_hard_partition)),
    method = "SM"
  )
  class <- max.col(consens_part$.Data)
  class <- match(class,unique(class))
  cl <- data.frame(part=factor(class), x=levels(se_sub$x))
  colData(se) <- merge(colData(se),cl,by="x")%>% DataFrame() %>% `row.names<-`(colnames(se_sub))
  metadata(se)$partition <- cl
  return(se)
}
