#' Quality control on cells
#'
#' @param sce SingleCellExperiment with \code{counts} and \code{ratio}
#' @param spike the character name of spike genes.
#' If missing, \code{spikePercent} will all be zero and
#' \code{filter_spike} will be false.
#' @param threshold A numeric scalar specifying the threshold
#' above which a gene is considered to be detected.
#' @param mad_sum A numeric scalar specifying exceed
#' how many median absolute deviations
#' from the median log10-counts a cell is considered to be filtered out.
#' Default is 5.
#' @param mad_detected A numeric scalar specifying exceed
#' how many median absolute deviations
#' from the median detected features a cell is considered to be filtered out.
#' Default is 5.
#' @param mad_spikegenes A numeric scalar specifying exceed
#' how many median absolute deviations
#' from the median spike genes expression percentage a cell
#' is considered to be filtered out. Default is 5.
#'
#' @return A DataFrame of QC statistics includes
#' \itemize{
#'   \item{sum} {the sum of counts of each cell}
#'   \item{detected} {the number of features above \code{threshold}}
#'   \item{spikePercent} {the percentage of counts assignes to spike genes}
#'   \item{filter_sum} {indicate whether log10-counts
#'    within \code{mad_sum} median absolute deviations
#' from the median log10-counts for the dataset}
#'   \item{filter_detected} {indicate whether features detected by this cell
#'    within \code{mad_detected} median absolute deviations
#' from the median detected features for the dataset}
#'   \item{filter_spike} {indicate whether percentage expressed by spike genes
#'    within \code{mad_spikegenes} median absolute deviations from the median
#'    spike genes expression percentage for the dataset}
#'   }
#'
#' @examples
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' cellQCmetrics <- cellQC(sce)
#' keep_cell <- (
#'   cellQCmetrics$filter_sum | # sufficient features (genes)
#'     # sufficient molecules counted
#'     cellQCmetrics$filter_detected |
#'     # sufficient features expressed compared to spike genes
#'     cellQCmetrics$filter_spike
#' )
#' sce <- sce[, keep_cell]
#'
#' # or manually setting threshold
#' cellQCmetrics <- cellQC(sce,
#'   spike = "Ercc",
#'   mad_detected = 4, mad_spikegenes = 4
#' )
#' keep_cell <- (
#'   cellQCmetrics$sum > 2.4 |
#'     cellQCmetrics$detected > 110
#' )
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment colData colData<- rowData
#' rowData<- assayNames assays assays<-
#' @importFrom S4Vectors DataFrame metadata metadata<-
#'
#' @export
cellQC <- function(sce, spike, threshold = 0,
                   mad_sum = 5, mad_detected = 3, mad_spikegenes = 5) {
  keep_feature <- rowSums(counts(sce)) > 0
  sce <- sce[keep_feature, ]
  sum <- log10(colSums(counts(sce)))
  filter_sum <- (sum > (median(sum) - mad_sum * mad(sum)))

  detected <- colSums(counts(sce) > threshold)
  filter_detected <- (detected > (median(detected) -
    mad_detected * mad(detected)))
  if (missing(spike)) {
    spikePercent <- rep(0, ncol(sce))
    filter_spike <- rep(TRUE, ncol(sce))
  } else {
    spike_gene <- grep(spike, row.names(sce))
    if (length(spike_gene) == 0) {
      message("No spike genes found")
    }
    spikePercent <- colSums(counts(sce)[spike_gene, ], na.rm = TRUE) * 100 /
      colSums(counts(sce)[-spike_gene, ], na.rm = TRUE)
    filter_spike <- (spikePercent > (median(spikePercent) -
      mad_spikegenes * mad(spikePercent)) &
      spikePercent < (median(spikePercent) +
        mad_spikegenes * mad(spikePercent)))
  }

  coldata <- cbind(
    colData(sce), sum, detected, spikePercent,
    filter_sum, filter_detected, filter_spike
  )
  coldata
}

#' Quality control on features
#'
#' @param sce SingleCellExperiment with \code{counts} and \code{ratio}
#' @param spike the character name of spike genes. The default is \code{Ercc}
#' @param detection_limit Numeric scalar providing the value above which
#' observations are deemed to be expressed.
#' @param threshold A numeric scalar specifying the threshold above which
#' percentage of cells expressed within each cell type. Default is 0.25
#' @param sd A numeric scalar specifying the cell type weighted
#' allelic ratio mean standard deviation threshold
#' above which are interested features with highly variation. Default is 0.03
#' @param pc pseudocount in the \code{preprocess} step
#'
#' @return A DataFrame of QC statistics includes
#' \itemize{
#'   \item{filter_celltype} {indicate whether genes expressed in more than
#' \code{threshold} cells for all cell types}
#'   \item{sd} {read counts standard deviation for each feature}
#'   \item{filter_sd} {indicate whether gene standard deviation
#' exceed \code{sd}}
#'   \item{filter_spike} {indicate no spike genes}
#'   }
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' featureQCmetric <- featureQC(sce)
#' keep_feature <- (featureQCmetric$filter_celltype &
#'   featureQCmetric$filter_sd &
#'   featureQCmetric$filter_spike)
#' sce <- sce[keep_feature, ]
#'
#' # or manually setting threshold
#' featureQCmetric <- featureQC(sce,
#'   spike = "Ercc",
#'   threshold = 0.25, sd = 0.03, pc = 2
#' )
#' keep_feature <- (featureQCmetric$filter_celltype &
#'   featureQCmetric$sd > 0.02)
#' @importFrom pbapply pbsapply
#' @importFrom scater nexprs
#'
#' @export
featureQC <- function(sce, spike, detection_limit = 1, threshold = 0.25, sd = 0.03, pc = 2) {
  check <- pbsapply(levels(sce$x), function(c) {
    poi <- which(sce$x == c)
    ct_threshold <- nexprs(counts(sce),
      byrow = TRUE, detection_limit = detection_limit,
      subset_col = poi
    ) >=
      length(poi) * threshold
    weighted_mean <- rowSums(assays(sce[, poi])[["ratio_pseudo"]] *
      (counts(sce[, poi]) + 2 * pc) /
      rowSums(counts(sce[, poi]) + 2 * pc))
    return(list(ct_threshold, weighted_mean))
  })
  ct_threshold <- do.call(cbind, check[seq(1, length(check), 2)])
  filter_celltype <- rowSums(ct_threshold) == nlevels(sce$x)

  weighted_mean <- do.call(cbind, check[seq(2, length(check), 2)])
  gsd <- rowSds(weighted_mean)
  filter_sd <- gsd > sd

  filter_spike <- rep(TRUE, nrow(sce))
  if (!missing(spike)) {
    filter_spike[grep(spike, row.names(sce))] <- FALSE
  }

  rowdata <- cbind(rowData(sce), filter_celltype, gsd, filter_sd, filter_spike)
  rowdata
}
