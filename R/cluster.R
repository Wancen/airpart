#' Gene clustering based on allelic ratio matrix with pseudo-count
#'
#' @param sce SingleCellExperiment containing assays \code{"ratio_pseudo"} and
#' colData factor \code{"x"}
#' @param method either \code{"GMM"} or \code{"hierarchical"}
#' @param plot logical, whether to make a PCA plot
#'
#' @references
#'
#' This function leverages Mclust from the mclust package, or hclust.
#'
#' For mclust see:
#' Luca Scrucca and Michael Fop and T. Brendan Murphy, Adrian E. Raftery
#' "mclust 5: clustering, classification and density
#' estimation using {G}aussian finite mixture models"
#' 2016. The R Journal. doi: 10.32614/RJ-2016-021
#'
#' @seealso \code{\link[mclust]{Mclust}}
#'
#' @return gene cluster IDs are stored in the rowData column \code{cluster}
#'
#' @importFrom mclust Mclust hc hcEII mclustBIC
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_minimal labs
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stats prcomp hclust cutree dist as.dist median mad
#' @importFrom grDevices rainbow
#' @importFrom graphics par
#'
#' @export
geneCluster <- function(sce, G = c(8, 12, 16, 20, 24),
                        method = "GMM", plot = TRUE, ...) {
  if (!"x" %in% names(colData(sce))) {
    stop('require a vector of annotated cell types "x" in colData')
  }
  if (!"ratio_pseudo" %in% assayNames(sce)) {
    stop('require an assay "ratio_pseudo"')
  }
  if (!is.factor(sce$x)) {
    sce$x <- as.factor(sce$x)
  }
  nct <- nlevels(sce$x)
  # PCA first
  pca <- prcomp(assays(sce)[["ratio_pseudo"]], rank. = 2 * nct) # use 2*nct
  ratio_pca <- as.matrix(pca$x)
  if (method == "GMM") {
    init <- list(hcPairs = mclust::hc(ratio_pca, modelName = "EII", use = "VARS"))
    d_clust <- mclust::Mclust(ratio_pca,
      G = G, modelNames = "EII",
      initialization = init, ...
    )
    nclust <- dim(d_clust$z)[2]
    my.clusters <- unname(d_clust$classification)
    cat("model-based optimal number of clusters:", nclust, "\n")
  }
  if (method == "hierarchical") {
    my.dist <- dist(ratio_pca, method = "manhattan")
    my.tree <- hclust(my.dist, method = "ward.D2")
    my.clusters <- unname(
      cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0)
    )
    nclust <- max(my.clusters)
  }
  if (plot) {
    dat <- as.data.frame(ratio_pca)
    dat$GeneCluster <- factor(my.clusters)
    cols <- unname(rainbow(nclust + 1)[-1])
    p <- ggplot(dat, aes(PC1, PC2, col = GeneCluster)) +
      geom_point() +
      scale_color_manual(values = cols) +
      theme_minimal()
    print(p)
  }
  rowdata <- cbind(rowData(sce), cluster = my.clusters)
  rowData(sce) <- rowdata
  metadata(sce)$geneCluster <- table(my.clusters)
  return(sce)
}
