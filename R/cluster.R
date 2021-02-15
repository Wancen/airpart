#' @import mclust
#' @import ggplot2
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export
#' @title Gene clustering based on allelic ratio matrix with pseudo count
#'
#' @description Clustering of genes
#'
#' @param ratio an allelic ratio matrix with pseudo count added
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
genecluster <- function(ratio, nct, G=c(4,8,12,16,20),
                        method="GMM", plot=TRUE,...) {
  # PCA first
  pca <- prcomp(ratio, rank. = 2*nct) #use 2*nct
  # vari<-summary(pca)$importance[2,1:2*nct]
  # sum(vari)
  ratio_pca <- as.matrix(pca$x)

  if (method=="GMM") {
    init <- list(hcPairs=hc(ratio_pca, modelName = "EII", use="VARS"))
    d_clust <- Mclust(ratio_pca, G=G, modelNames = "EII",
                      initialization=init, ...)
    # summary(d_clust)
    # plot(d_clust)
    nclust <- dim(d_clust$z)[2]
    my.clusters <- d_clust$classification
    cat("model-based optimal number of clusters:", nclust, "\n")
  }
  if (method=="hierarchical") {
    my.dist <- dist(ratio_pca, method="manhattan")
    my.tree <- hclust(my.dist, method="ward.D2")
    my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
    nclust <- max(my.clusters)
  }

  if (plot==TRUE) {
    p <- ggplot(data.frame(ratio_pca),
                aes(ratio_pca[,1], ratio_pca[,2],
                    col=as.factor(my.clusters))) +
      geom_point() +
      xlab("PC1") + ylab("PC2") +
      scale_color_manual(values=unname(palette.colors(nclust+1)[-1])) +
      theme_minimal() + labs(col = "Gene Cluster")
    print(p)
  }

  return(my.clusters)
}
