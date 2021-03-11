#' Make simulated data for airpart
#'
#' @param ngenecl number of genes per gene cluster
#' @param mu1 low count (typical of "noisy" ratio estimates)
#' @param mu2 high count
#' @param nct number of cell types
#' @param n number of cells per cell type
#' @param ngenecl number of genes per cluster
#' @param theta overdispersion parameter (higher is closer to binomial)
#'
#' @return SummarizedExperiment with the following elements as assays
#' \itemize{
#'   \item{ase.mat} {maternal allelic expression matrix}
#'   \item{ase.pat} {paternal allelic expression matrix}
#'   \item{true.ratio} {a matrix of the true probabilities
#' (allelic ratios) for the cell types}
#' }
#' Also \code{x} in the colData is a vector of annotated
#' cell types in the same order as cells in count matrix
#'
#' @importFrom emdbook rbetabinom
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
makeSimulatedData <- function(mu1, mu2, nct, n, ngenecl, theta){

  if ((nct %% 2) == 1) {
    stop("require even number of cell types for simulation setup")
  }

  ncl <- 3 # number of gene clusters
  ngene <- ncl * ngenecl # total number of genes
  nclcell <- nct * n * ngenecl # number elements within each gene cluster
  mean_total_count <- rep(rep(c(mu1, mu2),each=n/2), times=nct * ngene) # mean total count
  cts <- matrix(rpois(n * nct * ngene, mean_total_count), # total count matrix
                nrow = ngene, byrow = TRUE)

  # sets up the three gene clusters:
  # first cluster has pairs of cell types with same allelic ratio spanning from 1/10 to 9/10
  # second cluster has balanced allelic ratio
  # third cluster has pairs of cell types with same allelic ratio spanning from 7/10 to 9/10
  p.vec <- ( rep(c(
    seq(from=-4,to=4,length.out=nct/2),
    rep(0,nct/2),
    seq(from=2,to=4,length.out=nct/2)
  ), each=2) + 5 ) / 10
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

  true.ratio <- matrix(c(rep(p.vec[1:nct], ngenecl),
                         rep(p.vec[(nct+1):(2*nct)], ngenecl),
                         rep(p.vec[(2*nct+1):(3*nct)], ngenecl)),
                       ncol=nct,
                       byrow=TRUE)
  colnames(true.ratio) <- paste0("ct", 1:nct) # cell type names

  coldata <- data.frame(x=x)
  rowdata <- data.frame(true.ratio)
  assay.list <- list(ase=ase.mat, total=cts)
  SummarizedExperiment(assays=assay.list, colData=coldata, rowData=rowdata)
}
