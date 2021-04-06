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
#' @param p.vec the allelic ratio vector which follows gene cluster order.
#' (length is nct * ncl)
#'
#' @return SingleCellExperiment with the following elements as assays
#' \itemize{
#'   \item{a1} {allelic count matrix for the numerator/effect allele}
#'   \item{a2} {allelic count matrix for the denominator/non-effect allele}
#'   \item{true.ratio} {a matrix of the true probabilities
#' (allelic ratios) for the cell types}
#' }
#' Also \code{x} in the colData is a vector of annotated
#' cell types in the same order as cells in count matrix
#'
#' @examples
#' library(SummarizedExperiment)
#' sce <- makeSimulatedData()
#' assayNames(sce)
#' @importFrom emdbook rbetabinom
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats rnbinom
#'
#' @export
makeSimulatedData <- function(mu1 = 2, mu2 = 10, nct = 4, n = 30,
                              ngenecl = 50, theta = 20, ncl = 3,
                              p.vec = rep(c(0.2, 0.8, 0.5, 0.5, 0.7, 0.9),
                                each = 2
                              )) {
  if ((length(p.vec) / ncl) != nct) {
    stop("allelic ratio number is not matched with the product of
  number of cell types and number of gene cluster")
  }

  ngene <- ncl * ngenecl # total number of genes
  nclcell <- nct * n * ngenecl # number elements within each gene cluster
  # mean total count
  mean_total_count <- rep(rep(c(mu1, mu2), each = n / 2), times = nct * ngene)
  # total count matrix
  cts <- matrix(rnbinom(n * nct * ngene, mu = mean_total_count, size = 5),
    nrow = ngene, byrow = TRUE
  )

  p <- rep(p.vec, each = n * nct * ngene / length(p.vec))

  # allelic expression matrix for the effect allele
  a1 <- lapply(seq_len(ncl), function(m) {
    matrix(emdbook::rbetabinom(nclcell,
      prob = p[(nclcell * m - nclcell + 1):(nclcell * m)],
      size = cts[(m * ngenecl - ngenecl + 1):(m * ngenecl), ],
      theta = theta
    ), ncol = nct * n)
  })
  a1 <- do.call(rbind, a1)
  colnames(a1) <- paste0("cell", seq_len(nct * n))
  rownames(a1) <- paste0("gene", seq_len(ngene))
  a2 <- cts - a1

  x <- factor(rep(paste0("ct", seq_len(nct)), each = n)) # cell type vector

  true.ratio <- matrix(vapply(seq_len(ncl), function(m) {
    rep(p.vec[((m - 1) * nct + 1):(m * nct)], ngenecl)
  }, double(nct * ngenecl)),
  ncol = nct, byrow = TRUE
  )

  colnames(true.ratio) <- paste0("ct", seq_len(nct)) # cell type names
  coldata <- data.frame(x = factor(x, levels = unique(x)))
  rowdata <- data.frame(true.ratio)
  assay.list <- list(a1 = a1, a2 = a2)
  SingleCellExperiment(
    assays = assay.list,
    colData = coldata, rowData = rowdata
  )
}
