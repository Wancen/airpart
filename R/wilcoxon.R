#' Extension on Pairwise Mann Whitney Wilcoxon Test for partitioning
#'
#' Extends the Pairwise Mann Whitney Wilcoxon Test by combining
#' hierarchical clustering for partition.
#'
#' @param data A data frame containing the model response (\code{ratio}),
#' and grouping factors (\code{x})
#' @param threshold a vector with candidate thresholds for raw p-value
#' cut-off. For details please see vignette
#' @param p.adjust.method method for adjusting p-values
#' (see \code{\link[stats]{p.adjust}}). Can be abbreviated
#'
#' @return A vector grouping factor partition is returned
#' (TODO: recommend having same output as fusedlass()...x )
#'
#' @export
wilcox <- function(data,threshold=0.05,p.adjust.method="none",...) {
  nct <- length(levels(data$x))
  res  <-  pairwise.wilcox.test(data$ratio,data$x,p.adjust.method=p.adjust.method,...)
  # res <- pairwise.wilcox.test(dat$ratio,dat$x, p.adjust.method ="holm")
  adj <- as.data.frame(res$p.value)[lower.tri(res$p.value, diag = T)]
  adj <- ifelse(is.nan(adj),1,adj) # Wilcoxon output Nan if ratio of two cell types are exactly same
  b <- matrix(0, nct, nct)
  b[lower.tri(b, diag=FALSE)]=adj
  b2 <- b+t(b)
  diag(b2) <- 1
  bb <- ifelse(b2<threshold,1,0) # binarize p-value to be seen as dismilarity matrix
  clust <- hclust(as.dist(bb)) # hierarchical cluster on adjacency matrix
  my.clusters <- cutree(clust,h=0)
  return(my.clusters)
}

wilcox_adj <- function(data,threshold,p.adjust.method="none",...) {
  out <- list()
  obj <- sapply (1:length(threshold), function(j){
    fit <- wilcox(data,p.adjust.method=p.adjust.method,threshold=threshold[j],...)
    label <- tibble(type=factor(seq_along(1:nct)),par=fit)
    data2 <- data %>%
      left_join(label,by=c("x"="type"))
    data2 <- data2 %>%
      group_by(par) %>%
      mutate(grpmean=mean(ratio,na.rm = T))
    # loss function
    loss1 <- nrow(data) * log(sum((data2$ratio-data2$grpmean)^2,na.rm = T) /
                              nrow(data2))+length(unique(fit))*log(nrow(data2))
    out[["cl"]] <- fit
    out[["loss1"]] <- loss1
    return(out)
  })
  # TODO: can you fix the warning:
  # ...The `x` argument of `as_tibble.matrix()` must have unique column names...
  cl <- as_tibble(do.call(rbind, obj[seq(1,length(obj), by = 2)]))
  loss1 <- as_tibble(do.call(rbind, obj[seq(2,length(obj), by = 2)]))
  return(cl[which.min(unlist(loss1)),])
}

