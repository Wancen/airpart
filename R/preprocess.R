#' Preprocess the SingleCellExperiment
#'
#' @param sce SingleCellExperiment with \code{a1} (effect allele)
#' and \code{a2} (non-effect allele). The allelic ratio will be
#' calculated as a1 / (a1 + a2).
#' @param pc pseudocount for calculating the smoothed ratio
#'
#' @return SingleCellExperiment with total count,
#' allelic ratio = a1/(a1 + a2), and pseud-ocount-smoothed ratio
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
  if (any(assays(sce)[["a1"]] != round(assays(sce)[["a1"]]))) {
    warning("Functions estDisp, fusedLasso expect integers")
  }
  assays(sce)[["counts"]] <- assays(sce)[["a1"]] +
    assays(sce)[["a2"]]
  assays(sce)[["ratio"]] <- assays(sce)[["a1"]] /
    assays(sce)[["counts"]]
  assays(sce)[["ratio_pseudo"]] <- (assays(sce)[["a1"]] + pc) /
    (assays(sce)[["counts"]] + 2 * pc)
  metadata(sce) <- list(airpartVersion = packageVersion("airpart"))
  sce
}
