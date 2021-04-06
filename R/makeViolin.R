#' Smoothed allelic ratio in violin plots
#'
#' @param sce SingleCellExperiment
#'
#' @return a ggplot2 object
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' makeViolin(sce_sub)
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_violin
#' scale_fill_brewer theme labs
#' @importFrom dplyr n
#'
#' @export
makeViolin <- function(sce) {
  cl_ratio <- as.vector(unlist(assays(sce)[["ratio_pseudo"]]))
  cl_total <- as.vector(unlist(counts(sce)))
  dat <- data.frame(
    ratio = cl_ratio,
    x = factor(rep(sce$x, each = length(sce))),
    cts = cl_total,
    part = factor(rep(sce$part, each = length(sce)))
  )
  dat <- dat[!is.nan(dat$ratio), ]
  # sample size
  sample_size <- dat %>%
    group_by(.data$x) %>%
    summarise(num = n())
  dat <- dat %>%
    left_join(sample_size)
  dat$myaxis <- paste0(dat$x, "\n", "n=", dat$num)
  p <- ggplot(dat, aes(x = .data$myaxis, y = .data$ratio, fill = .data$part)) +
    geom_violin(width = 1.2, color = "#A4A4A4") +
    geom_boxplot(width = 0.1, color = "black") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "BuPu") +
    labs(x = "cell type", y = "allelic ratio")
  p
}
