#' Consensus Partitions
#'
#' Derive consensus partitions of an ensemble fused lasso partitions.
#'
#' @param sce SingleCellExperiment
#'
#' @return A matrix grouping factor partition is replaced in metadata.
#' Consensus Partation also stored in colData\code{"part"}.
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G=1:4)
#' f <- ratio ~ p(x, pen = "gflasso") # formula for the GFL
#' sce_sub <- fusedLasso(sce,formula=f,model="binomial",
#'                       genecluster=1,ncores=4, niter=2)
#' sce_sub <- consensusPart(sce_sub)
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
