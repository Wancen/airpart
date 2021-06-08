#' Allelic ratio summary
#'
#' Oroduce allelic ratio summaries for each gene cluster
#'
#' @param sce SingleCellExperiment
#' @param genecluster an optional vector of gene cluster IDs.
#' if nothing is given, all cluster's summaries will be
#' calculated
#'
#' @return a list of gene cluster summary tables containing:
#' \itemize{
#'   \item{weighted.mean} {weighted mean of allelic ratio for the cell types}
#'   \item{mean} {mean allelic ratio for the cell types}
#'   \item{var} {variance of allelic ratio for the cell types}
#' }
#'
#' @examples
#'
#' library(S4Vectors)
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' summary <- summaryAllelicRatio(sce, genecluster = c(1, 3))
#' summary
#' @importFrom dplyr summarise group_by %>%
#' @importFrom stats var weighted.mean
#'
#' @export
summaryAllelicRatio <- function(sce, genecluster) {
  if (missing(genecluster)) {
    genecluster <- seq_len(max(rowData(sce)$cluster))
  }
  res <- lapply(genecluster, function(i) {
    sce_sub <- sce[rowData(sce)$cluster == i, ]
    cl_ratio <- as.vector(unlist(assays(sce_sub)[["ratio"]]))
    cl_total <- as.vector(unlist(assays(sce_sub)[["counts"]]))
    dat <- data.frame(
      ratio = cl_ratio,
      x = factor(rep(sce_sub$x, each = length(sce_sub))),
      cts = cl_total
    )
    summary <- dat %>%
      group_by(.data$x) %>%
      summarise(
        weighted.mean = weighted.mean(.data$ratio, .data$cts,
          na.rm = TRUE
        ),
        mean = mean(.data$ratio, na.rm = TRUE),
        var = var(.data$ratio, na.rm = TRUE)
      ) %>%
      as.data.frame()
    summary
  })
  names(res) <- paste(
    "gene cluster", genecluster, "with",
    metadata(sce)$geneCluster[genecluster], "genes"
  )
  return(res)
}
