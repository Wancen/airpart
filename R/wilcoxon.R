#' Extension on Pairwise Mann Whitney Wilcoxon Test for partitioning
#'
#' Extends the Pairwise Mann Whitney Wilcoxon Test by combining
#' hierarchical clustering for partition.
#'
#' @param sce A SingleCellExperiment containing assays (\code{"ratio"},
#' \code{"counts"}) and colData \code{"x"}
#' @param genecluster which gene cluster result want to be returned.
#' Usually identified interesting gene cluster pattern by \code{\link{summaryAllelicRatio}}
#' @param threshold a vector with candidate thresholds for raw p-value
#' cut-off. For details please see vignette
#' @param p.adjust.method method for adjusting p-values
#' (see \code{\link[stats]{p.adjust}}). Can be abbreviated
#' @param ... additional arguments to pass to \code{\link[stats]{wilcox.test}}.
#'
#' @return A vector grouping factor partition is returned
#'
#' @importFrom dplyr left_join
#' @importFrom plyr mutate
#' @importFrom stats pairwise.wilcox.test
#'
#' @export
wilcoxExt <- function(sce, genecluster, threshold, p.adjust.method = "none", ...) {
  # construct data frame
  sce_sub <- sce[rowData(sce)$cluster == genecluster, ]
  cl_ratio <- as.vector(unlist(assays(sce_sub)[["ratio"]]))
  cl_total <- as.vector(unlist(counts(sce_sub)))
  dat <- data.frame(
    ratio = cl_ratio,
    x = factor(rep(sce_sub$x, each = length(sce_sub))),
    cts = cl_total
  )
  nct <- nlevels(sce$x)

  out <- list()
  obj <- sapply(1:length(threshold), function(j) {
    fit <- wilcoxInt(dat, p.adjust.method = p.adjust.method, threshold = threshold[j],...)
    label <- data.frame(type = factor(seq_along(1:nct)), par = fit)
    dat2 <- dat %>%
      left_join(label, by = c("x" = "type"))
    dat2 <- dat2 %>%
      group_by(par) %>%
      mutate(grpmean = mean(ratio, na.rm = TRUE))
    # loss function
    loss1 <- nrow(dat) * log(sum((dat2$ratio - dat2$grpmean)^2, na.rm = TRUE) /
      nrow(dat2)) + length(unique(fit)) * log(nrow(dat2))
    out[["cl"]] <- fit
    out[["loss1"]] <- loss1
    return(out)
  })

  cl <- do.call(rbind, obj[seq(1, length(obj), by = 2)])
  loss1 <- do.call(rbind, obj[seq(2, length(obj), by = 2)])
  partition <- data.frame(part = factor(cl[which.min(loss1), ]), x = levels(sce_sub$x))
  coldata<-DataFrame(rowname=colnames(sce_sub), colData(sce_sub))
  coldata <- merge(coldata, partition, by = "x",sort=F) %>%
    DataFrame()
  rownames(coldata)<-coldata$rowname
  colData(sce_sub)<-coldata
  metadata(sce_sub)$partition <- partition
  return(sce_sub)
}

wilcoxInt <- function(data, threshold = 0.05, p.adjust.method = "none", ...) {
  nct <- length(levels(data$x))
  res <- pairwise.wilcox.test(data$ratio, data$x, p.adjust.method = p.adjust.method, ...)
  adj <- as.data.frame(res$p.value)[lower.tri(res$p.value, diag = TRUE)]
  adj <- ifelse(is.nan(adj), 1, adj) # Wilcoxon output Nan if ratio of two cell types are exactly same
  b <- matrix(0, nct, nct)
  b[lower.tri(b, diag = FALSE)] <- adj
  b2 <- b + t(b)
  diag(b2) <- 1
  bb <- ifelse(b2 < threshold, 1, 0) # binarize p-value to be seen as dismilarity matrix
  clust <- hclust(as.dist(bb)) # hierarchical cluster on adjacency matrix
  my.clusters <- cutree(clust, h = 0)
  return(my.clusters)
}
