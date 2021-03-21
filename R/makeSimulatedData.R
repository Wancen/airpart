#' Make simulated data for airpart
#'
#' @param ngenecl number of genes per gene cluster
#' @param mu1 low count (typical of "noisy" ratio estimates)
#' @param mu2 high count
#' @param nct number of cell types
#' @param n number of cells per cell type
#' @param ngenecl number of genes per cluster
#' @param theta overdispersion parameter (higher is closer to binomial)
#' @param ncl number of gene cluster
#' @param p.vec the allelic ratio vector which follows gene cluster order. (length is nct * ncl)
#'
#' @return SingleCellExperiment with the following elements as assays
#' \itemize{
#'   \item{ase.mat} {maternal allelic expression matrix}
#'   \item{ase.pat} {paternal allelic expression matrix}
#'   \item{true.ratio} {a matrix of the true probabilities
#' (allelic ratios) for the cell types}
#' }
#' Also \code{x} in the colData is a vector of annotated
#' cell types in the same order as cells in count matrix
#'
#' @examples
#'
#' library(SummarizedExperiment)
#' sce <- makeSimulatedData()
#' assayNames(sce)
#'
#' @importFrom emdbook rbetabinom
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats rnbinom
#'
#' @export
makeSimulatedData <- function(mu1 = 2, mu2 = 10, nct = 4, n = 30,
                              ngenecl = 50, theta = 20, ncl = 3,
                              p.vec = rep(c(0.2,0.8,0.5,0.5,0.7,0.9),each = 2)){

  if ((length(p.vec) / ncl) != nct) {
    stop("allelic ratio number is not matched with the product of number of cell types and number of gene cluster")
  }

  ngene <- ncl * ngenecl # total number of genes
  nclcell <- nct * n * ngenecl # number elements within each gene cluster
  mean_total_count <- rep(rep(c(mu1, mu2),each=n/2), times=nct * ngene) # mean total count
  cts <- matrix(rnbinom(n * nct * ngene, mu=mean_total_count,size=2), # total count matrix
                nrow = ngene, byrow = TRUE)

  p <- rep(p.vec, each=n*nct*ngene/length(p.vec))

  # maternal allelic expression matrix
  ase.mat <- lapply(1:ncl,function(m) {
    matrix(emdbook::rbetabinom(nclcell, prob=p[(nclcell*m-nclcell+1):(nclcell*m)],
                               size=cts[(m*ngenecl-ngenecl+1):(m*ngenecl),],
                               theta=theta),ncol = nct*n)})
  ase.mat <- do.call(rbind,ase.mat)
  colnames(ase.mat) <- paste0("cell",1:(nct*n))
  rownames(ase.mat) <- paste0("gene",1:ngene)

  ase.pat<-cts-ase.mat # paternal allelic expression matrix

  x <- factor(rep(paste0("ct", 1:nct),each=n)) # cell type vector

  true.ratio<-matrix(sapply(1:ncl, function(m){
    rep(p.vec[((m-1)*nct+1):(m*nct)], ngenecl)
  }),
  ncol=nct,byrow=TRUE)

  colnames(true.ratio) <- paste0("ct", 1:nct) # cell type names
  coldata <- data.frame(x=factor(x,levels = unique(x)))
  rowdata <- data.frame(true.ratio)
  assay.list <- list(ase.mat=ase.mat, ase.pat=ase.pat)
  SingleCellExperiment(assays=assay.list, colData=coldata, rowData=rowdata)
}
