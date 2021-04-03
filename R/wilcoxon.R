#' Extension on Pairwise Mann Whitney Wilcoxon Test for partitioning
#'
#' Extends the Pairwise Mann Whitney Wilcoxon Test by combining
#' hierarchical clustering for partition.
#'
#' @param sce A SingleCellExperiment containing assays (\code{"ratio"},
#' \code{"counts"}) and colData \code{"x"}
#' @param genecluster which gene cluster result want to be returned.
#' Usually identified interesting gene cluster pattern by
#' \code{\link{summaryAllelicRatio}}
#' @param threshold a vector with candidate thresholds for raw p-value
#' cut-off. Default is 10^seq(from=-2,to=-0.4,by=0.2).
#' For details please see vignette
#' @param p.adjust.method method for adjusting p-values
#' (see \code{\link[stats]{p.adjust}}). Can be abbreviated
#' @param ... additional arguments to pass to \code{\link[stats]{wilcox.test}}.
#' @param adj.matrix an adjacency matrix with 1 indicates cell states
#' allowed to be grouped together, 0 otherwise.
#'
#' @return A matrix grouping factor partition and
#' the significant cut-off threshold
#' are returned in metadata \code{"partition"} and \code{"threshold"}.
#' Partation also stored in colData\code{"part"}.
#'
#' @examples
#' library(S4Vectors)
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = seq_len(4))
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' metadata(sce_sub)$partition
#' metadata(sce_sub)$threshold
#'
#' # Suppose we have 4 cell states, if we don't want cell state 1
#' # to be grouped together with other cell states
#' adj.matrix <- 1 - diag(4)
#' colnames(adj.matrix) <- rownames(adj.matrix) <- levels(sce$x)
#' adj.matrix[1, c(2, 3, 4)] <- 0
#' adj.matrix[c(2, 3, 4), 1] <- 0
#' thrs <- 10^seq(from = -2, to = -0.4, by = 0.1)
#' sce_sub <- wilcoxExt(sce,
#'   genecluster = 1, threshold = thrs,
#'   adj.matrix = adj.matrix
#' )
#' metadata(sce_sub)$partition
#' @importFrom dplyr left_join
#' @importFrom plyr mutate
#' @importFrom stats pairwise.wilcox.test
#'
#' @export
wilcoxExt <- function(sce, genecluster, threshold,
                      p.adjust.method = "none", adj.matrix, ...) {
  if (missing(threshold)) {
    threshold <- 10^seq(from = -2, to = -0.4, by = 0.2)
  }
  if (missing(genecluster)) {
    stop("No gene cluster number")
  }

  stopifnot(c("ratio", "counts") %in% assayNames(sce))
  stopifnot("x" %in% names(colData(sce)))
  stopifnot("cluster" %in% names(rowData(sce)))

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

  if (missing(adj.matrix)) {
    adj.matrix <- matrix(1, nct, nct)
  }

  obj <- lapply(seq_len(length(threshold)), function(j) {
    fit <- wilcoxInt(dat,
      p.adjust.method = p.adjust.method,
      threshold = threshold[j], adj.matrix = adj.matrix, ...
    )
    label <- data.frame(type = factor(levels(sce$x)), par = factor(fit))
    dat2 <- dat %>%
      left_join(label, by = c("x" = "type"))
    dat2 <- dat2 %>%
      group_by(.data$par) %>%
      dplyr::mutate(grpmean = mean(.data$ratio, na.rm = TRUE))
    # loss function
    loss1 <- nrow(dat) * log(sum((dat2$ratio - dat2$grpmean)^2, na.rm = TRUE) /
      nrow(dat2)) + length(unique(fit)) * log(nrow(dat2))
    return(list(cl = fit, loss1 = loss1))
  })

  cl <- do.call(rbind, lapply(obj, `[[`, 1))
  loss1 <- do.call(rbind, lapply(obj, `[[`, 2))
  partition <- data.frame(part = factor(cl[which.min(loss1), ]),
                          x = levels(sce_sub$x))
  cd <- colData(sce_sub)
  cd2 <- cd[, !names(cd) %in% c("part", "rowname")] %>%
    as.data.frame() %>%
    setNames(names(cd)[!names(cd) %in% c("part", "rowname")])
  coldata <- DataFrame(rowname = colnames(sce_sub), cd2)
  coldata <- merge(coldata, partition, by = "x", sort = FALSE) %>%
    DataFrame()
  rownames(coldata) <- coldata$rowname
  colData(sce_sub) <- coldata
  metadata(sce_sub)$partition <- partition
  metadata(sce_sub)$threshold <- threshold[which.min(loss1)]
  return(sce_sub)
}

# TODO can we remove this code below then? 

# not exported
wilcoxInt <- function(data, threshold = 0.05,
                      p.adjust.method = "none", adj.matrix, ...) {
  nct <- length(levels(data$x))
  res <- pairwise.wilcox.test(data$ratio, data$x,
                              p.adjust.method = p.adjust.method, ...)
  adj <- as.data.frame(res$p.value)[lower.tri(res$p.value, diag = TRUE)]
  # Wilcoxon output Nan if ratio of two cell types are exactly same
  adj <- ifelse(is.nan(adj), 1, adj) 
  b <- matrix(0, nct, nct)
  b[lower.tri(b, diag = FALSE)] <- adj
  b2 <- b + t(b)
  diag(b2) <- 1
  b2[which(adj.matrix == 0)] <- 0
  # binarize p-value to be seen as dismilarity matrix
  bb <- ifelse(b2 < threshold, 1, 0)
  # hierarchical cluster on adjacency matrix
  clust <- hclust(as.dist(bb)) 
  my.clusters <- cutree(clust, h = 0)
  return(my.clusters)
}
