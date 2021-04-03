#' Preprocess the SingleCellExperiment
#'
#' @param sce SingleCellExperiment with \code{ase.mat} and \code{ase.pat}
#' @param pc pseudocount for calculating the smoothed ratio
#'
#' @return SingleCellExperiment with total count,
#' ratio and pseudocount-smoothed ratio, where the
#' ratio provides the maternal allele count over total
#'
#' @examples
#' library(SummarizedExperiment)
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' assayNames(sce)
#' @importFrom utils packageVersion
#'
#' @export
preprocess <- function(sce, pc = 2) {
  assays(sce)[["counts"]] <- assays(sce)[["ase.mat"]] +
    assays(sce)[["ase.pat"]]
  assays(sce)[["ratio"]] <- assays(sce)[["ase.mat"]] /
    assays(sce)[["counts"]]
  assays(sce)[["ratio_pseudo"]] <- (assays(sce)[["ase.mat"]] + pc) /
    (assays(sce)[["counts"]] + 2 * pc)
  metadata(sce) <- list(airpartVersion = packageVersion("airpart"))
  sce
}
