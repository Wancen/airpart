#' Posterior mean allelic ratio estimates in violin plots
#'
#' @param sce SingleCellExperiment
#' @param xlab the x axis name.
#' @param ylim the y axis range
#'
#' @return a ggplot2 object, \code{n} represents number of cells in that cell type.
#'
#' @examples
#'
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' sce_sub <- allelicRatio(sce_sub)
#' makeViolin(sce_sub)
#' @export
makeViolin <- function(sce, xlab = "cell type", ylim = c(0,1)) {
  ylim1<- ylim[1]
  ylim2<- ylim[2]
  ar <- rowData(sce)[, c(grep("ar_", colnames(rowData(sce)), value = TRUE))] %>%
    `colnames<-`(levels(sce$x))
  dat <- data.frame(
    ratio = as.vector(unlist(ar)),
    x = factor(rep(levels(sce$x), each = length(sce))),
    part = factor(rep(metadata(sce)$partition$part, each = length(sce)))
  )
  sample_size <- data.frame(table(sce$x)) %>% `colnames<-`(c("x", "n"))
  dat <- dat %>%
    left_join(sample_size)
  dat$myaxis <- paste0(dat$x, "\n", "n=", dat$n) %>% as.factor()
  dat$myaxis <- factor(dat$myaxis, levels = unique(paste0(dat$x, "\n", "n=", dat$n)))
  p <- ggplot(dat, aes(x = .data$myaxis, y = .data$ratio, fill = .data$part)) +
    geom_violin(color = "#A4A4A4", size = 1.2, alpha = .7) +
    geom_boxplot(width = 0.1, color = "black") +
    theme_minimal() + ylim(ylim1,ylim2) + 
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = "transparent")
    ) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = xlab, y = "allelic ratio")
  p
}
