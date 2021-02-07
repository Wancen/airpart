#' @import mclust
#' @import ggplot2
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export
#' @title Gene clustering based on allelic ratio matrix with pseudo count
#'
#' @description cluster on genes
#'
#' @param ratio an allelic ratio matrix with pseudo count added in
#' @param nct number of cell types or states
#' @param method either \code{"GMM"} or \code{"hierarchical"}
genecluster<-function(ratio,nct,G=c(4,8,12,16,20),method="GMM",plot=TRUE,...){
  #PCA first
  pca<-prcomp(ratio,rank. = 2*nct) #use 2*nct
  # vari<-summary(pca)$importance[2,1:2*nct]
  # sum(vari)
  ratio_pca<-as.matrix(pca$x)

  if(method=="GMM"){
    d_clust<-Mclust(ratio_pca,G=G,modelNames = "EII",initialization = list(hcPairs = hc(ratio_pca,modelName = "EII", use = "VARS")),...)
    # summary(d_clust)
    # plot(d_clust)
    m.best <- dim(d_clust$z)[2]
    cat("model-based optimal number of clusters:", m.best, "\n")
    if(plot==TRUE){
      p<-ggplot(data.frame(ratio_pca),aes(ratio_pca[,1],ratio_pca[,2],col=as.factor(d_clust$classification))) +
        geom_point() +xlab("PC1")+ylab("PC2")+scale_color_manual(values = rainbow(m.best))+theme_minimal()+labs(col = "Gene Cluster")
      print(p)
    }
    # mod<-MclustDR(d_clust)
    # plot(mod)
    return(d_clust$classification)
  }
  if(method=="hierarchical"){
    my.dist <- dist(ratio_pca,method = "manhattan")
    my.tree <- hclust(my.dist, method="ward.D2")
    my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
    if(plot==TRUE){
      p<-ggplot(data.frame(ratio_pca),aes(ratio_pca[,1],ratio_pca[,2],col=as.factor(my.clusters))) +
        geom_point() +xlab("PC1")+ylab("PC2")+scale_color_manual(values = rainbow(m.best))+theme_minimal()+labs(col = "Gene Cluster")
      print(p)
    }
    return(my.clusters)
  }
}
