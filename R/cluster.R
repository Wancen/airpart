#' Gene clustering based on allelic ratio matrix with pseudo-count
#'
#' @param ratio an allelic ratio matrix with pseudo-count added
#' @param nct number of cell types or states
#' @param method either \code{"GMM"} or \code{"hierarchical"}
#' @param plot logical, whether to make a plot
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
#' @return a vector of the gene cluster IDs
#'
#' @importFrom mclust Mclust hc hcEII mclustBIC
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_minimal labs
#' @importFrom dynamicTreeCut cutreeDynamic
#'
#' @export
genecluster <- function(ratio, nct, G=c(4,8,12,16,20),
                        method="GMM", plot=TRUE,...) {
  # PCA first
  pca <- prcomp(ratio, rank. = 2*nct) # use 2*nct
  ratio_pca <- as.matrix(pca$x)
  if (method=="GMM") {
    init <- list(hcPairs=mclust::hc(ratio_pca, modelName="EII", use="VARS"))
    d_clust <- mclust::Mclust(ratio_pca, G=G, modelNames="EII",
                              initialization=init, ...)
    nclust <- dim(d_clust$z)[2]
    my.clusters <- d_clust$classification
    cat("model-based optimal number of clusters:", nclust, "\n")
  }
  if (method=="hierarchical") {
    my.dist <- dist(ratio_pca, method="manhattan")
    my.tree <- hclust(my.dist, method="ward.D2")
    my.clusters <- unname(
      cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0)
    )
    nclust <- max(my.clusters)
  }
  if (plot==TRUE) {
    dat <- as.data.frame(ratio_pca)
    dat$GeneCluster <- factor(my.clusters)
    cols <- unname(rainbow(nclust+1)[-1])
    p <- ggplot(dat, aes(PC1, PC2, col=GeneCluster)) +
      geom_point() +
      scale_color_manual(values=cols) +
      theme_minimal()
    print(p)
  }
  return(my.clusters)
}
