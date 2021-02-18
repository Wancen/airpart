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
#' @return a list with the following elements:
#' \itemize{
#'   \item{ase.mat} {maternal allelic expression matrix}
#'   \item{ase.pat} {paternal allelic expression matrix}
#'   \item{x} {a vector of annotated cell types in the same order as cells in count matrix}
#'   \item{true.ratio} {a matrix of the true probabilities (allelic ratios) for the cell types}
#' }
#'
#' @importFrom emdbook rbetabinom
#'
#' @export
makeSimulatedData <- function(mu1, mu2, nct, n, ngenecl, theta){

  if ((nct %% 2) == 1) {
    stop("require even number of cell types for simulation setup")
  }

  ncl <- 3 # number of gene clusters
  ngene <- ncl*ngenecl # total number of genes
  nclcell <- nct*n*ngenecl # number elements within each gene cluster
  mean_total_count <- rep(rep(c(mu1, mu2),each=n/2), times=nct*ngene) # mean total count
  cts <- matrix(rpois(n * nct*ngene, mean_total_count), # total count matrix
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
  ase.mat<-lapply(1:ncl,function(m) {
    matrix(rbetabinom(nclcell, prob=p[(nclcell*m-nclcell+1):(nclcell*m)],
                      size=cts[(m*ngenecl-ngenecl+1):(m*ngenecl),],
                      theta=theta),ncol = nct*n)})
  ase.mat <- do.call(rbind,ase.mat)
  colnames(ase.mat) <- paste0("cell",1:(nct*n))
  rownames(ase.mat) <- paste0("gene",1:ngene)

  ase.pat<-cts-ase.mat # paternal allelic expression matrix

  x <- factor(rep(1:nct,each=n)) # cell type vector

  list(ase.mat = ase.mat, ase.pat = ase.pat, x = x,
       true.ratio = matrix(p.vec,ncol=nct,byrow=TRUE,dimnames=list(1:3,1:nct)))
}
