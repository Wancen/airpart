#' Estimate overdispersion parameter of a betabinomial
#'
#' @param se SummarizedExperiment with \code{ase} and \code{total}
#' @param ct the name of the cell type to be estimated. Default is the highest one.
#'
#' @importFrom apeglm apeglm bbEstDisp
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_minimal labs

estDisp<-function(se, ct){
  if (missing(ct)) {
    cell_sum<-colSums(assays(se)[["total"]]) %>% by(se$x,FUN=mean)
    ct=names(which.max(cell_sum))
  }

  se_sub<-se[,se$x==ct]
  x <- matrix(1, ncol=1, nrow=dim(se_sub)[2])
  theta.hat <- 100
  param <- cbind(theta.hat, (assays(se_sub)[["total"]]+4))
  for (i in 1:5) {
    fit.mle <- apeglm(Y=(assays(se_sub)[["ase"]]+2), x=x, log.lik=NULL, param=param, no.shrink=TRUE, log.link=FALSE, method="betabinC")
    theta.hat <- bbEstDisp(success=(assays(se_sub)[["ase"]]+2), size=(assays(se_sub)[["total"]]+4), x=x, beta=fit.mle$map, minDisp=.01, maxDisp=5000)
  }
  gene_mean<-rowMeans(assays(se_sub)[["total"]])
  est<-data.frame(mean=gene_mean,theta=theta.hat)
  p<-ggplot(est[which(est$mean>2),],aes(mean, theta)) +
    geom_point() +geom_smooth()+ylim(0,75)+
    theme_minimal()+
    labs(x="gene mean",y="theta")
  print(p)
}




