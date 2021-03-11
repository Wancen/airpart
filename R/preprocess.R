#' Preprocess the SummarizedExperiment
#'
#' @param se SummarizedExperiment with \code{ase} and \code{total}
#' @param pc pseudocount for calculating the
#' smoothed ratio
#'
#' @return SummarizedExperiment with total count,
#' ratio and pseudocount-smoothed ratio, where the
#' ratio provides the maternal allele count over total
#'
#' @export
preprocess <- function(se, pc=2) {
  # assays(se)[["total"]] <- assays(se)[["ase.mat"]] + assays(se)[["ase.pat"]]
  assays(se)[["ratio"]] <- assays(se)[["ase"]] / assays(se)[["total"]]
  assays(se)[["ratio_pseudo"]] <- (assays(se)[["ase"]] + pc) /
    (assays(se)[["total"]] + 2*pc)
  metadata(se) <- list(airpartVersion=packageVersion("airpart"))
  se
}
