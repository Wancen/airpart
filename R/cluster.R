#' Gene clustering based on allelic ratio matrix with pseudo-count
#'
#' @param sce SingleCellExperiment containing assays \code{"ratio_pseudo"} and
#' colData factor \code{"x"}
#' @param method the method to do gene clustering. The default is the Gaussian Mixture Modeling which
#' is likely to be more accurate. \code{"hierarchical"} represents
#' automatic hierarchical clustering which is faster to compute.
#' @param plot logical, whether to make a PCA plot
#' @param G An integer vector specifying the numbers of clusters for which the BIC is to be calculated.
#'  The default is G=c(8, 12, 16, 20, 24).
#' @param ... Catches unused arguments in indirect or list calls via do.call
#' as described in \code{\link[mclust]{Mclust}}
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
#' and a table of gene cluster is returned in metadata \code{geneCluster}
#'
#' @examples
#' 
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' 
#' @importFrom mclust Mclust hc hcEII mclustBIC
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_minimal labs
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stats prcomp hclust cutree dist as.dist median mad
#' @importFrom grDevices rainbow
#' @importFrom graphics par
#' @importFrom rlang .data
#'
#' @export
geneCluster <- function(sce, G, method = c("GMM", "hierarchical"),
                        plot = TRUE, ...) {
  method <- match.arg(method, c("GMM", "hierarchical"))[1]
  if (missing(method)) {
    G <- c(8, 12, 16, 20, 24)
  }
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
    fit <- gmmCluster(ratio_pca, G, ...)
  } else if (method == "hierarchical") {
    fit <- hierCluster(ratio_pca)
  }
  nclust <- fit$clust
  my_clusters <- fit$my_clusters

  if (plot) {
    dat <- as.data.frame(ratio_pca)
    dat$GeneCluster <- factor(my_clusters)
    cols <- unname(rainbow(nclust + 1)[-1])
    p <- ggplot(dat, aes(.data$PC1, .data$PC2, col = .data$GeneCluster)) +
      geom_point() +
      scale_color_manual(values = cols) +
      theme_minimal()
    print(p)
  }
  rowdata <- cbind(rowData(sce), cluster = my_clusters)
  rowData(sce) <- rowdata
  metadata(sce)$geneCluster <- table(my_clusters)
  return(sce)
}

gmmCluster <- function(ratio_pca, G, ...) {
  init <- list(hcPairs = mclust::hc(ratio_pca, modelName = "EII", use = "VARS"))
  d_clust <- mclust::Mclust(ratio_pca,
    G = G, modelNames = "EII",
    initialization = init, ...
  )
  nclust <- dim(d_clust$z)[2]
  my_clusters <- unname(d_clust$classification)
  cat("model-based optimal number of clusters:", nclust, "\n")
  list(nclust = nclust, my_clusters = my_clusters)
}

hierCluster <- function(ratio_pca) {
  my_dist <- dist(ratio_pca, method = "manhattan")
  my_tree <- hclust(my_dist, method = "ward.D2")
  my_clusters <- unname(
    cutreeDynamic(my_tree, distM = as.matrix(my_dist), verbose = 0)
  )
  nclust <- max(my_clusters)
  list(nclust = nclust, my_clusters = my_clusters)
}
