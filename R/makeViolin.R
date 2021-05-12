#' Posterior mean allelic ratio estimates in violin plots
#'
#' @param sce SingleCellExperiment
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
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_violin
#' scale_fill_brewer theme labs
#'
#' @export
makeViolin <- function(sce) {
    ar <- rowData(sce)[, c(grep("ar", colnames(rowData(sce)), value = TRUE))] %>%
        `colnames<-`(levels(sce$x))
    dat <- data.frame(
        ratio = as.vector(unlist(ar)),
        x = factor(rep(levels(sce$x), each = length(sce))),
        part = factor(rep(metadata(sce)$partition$part, each = length(sce)))
    )
    sample_size <- data.frame(table(sce$x)) %>% `colnames<-`(c("x", "n"))
    dat <- dat %>%
        left_join(sample_size)
    dat$myaxis <- paste0(dat$x, "\n", "n=", dat$n)
    p <- ggplot(dat, aes(x = .data$myaxis, y = .data$ratio, fill = .data$part)) +
        geom_violin(width = 1.2, color = "#A4A4A4") +
        geom_boxplot(width = 0.1, color = "black") +
        theme_minimal() +
        theme(legend.position = "none") +
        scale_fill_brewer(palette = "Set2") +
        labs(x = "cell type", y = "allelic ratio")
    p
}
