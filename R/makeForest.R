#' Plot allelic ratio result as forest
#'
#' Draw a forest plot to visualized cell type specific
#' allelic ratio estimator and confidence interval.
#' It is based on the \pkg{forestplot}-package`s \code{forestplot} function.
#'
#' @param sce A SingleCellExperiment containing colData
#' allelic ratio estimator in the third column
#' and last two column is the confidence interval.
#' @param genepoi the gene position index or gene name vector that want to be plotted.
#' Ordered by increased cell type svalue.
#' Default is the top 40 genes that has minimum svalue in any cell type
#' or all genes if number of genes smaller than 40.
#' @param ctpoi the cell type position index that want to be plotted.
#' @param showtext indicate whether show the svalue information along the forestplot.
#' @param xticks argument as described in \code{\link[forestplot]{forestplot}}
#' @param boxsize Override the default box size based on precision
#' @param xlab x-axis label. Default is "Allelic Ratio"
#' @param col Set the colors for all the elements. See
#' \code{\link[forestplot]{fpColors}} for details
#' @param grid If you want a discrete gray dashed grid
#' at the level of the ticks
#' you can set this parameter to TRUE.
#' If you set the parameter to a vector of values lines will be drawn at the
#' corresponding positions.
#' If you want to specify the \code{\link[grid]{gpar}} of the lines
#' then either directly pass a \code{\link[grid]{gpar}} object
#' or set the gp attribute e.g.
#' \code{attr(line_vector, "gp") <- }
#' \code{\link[grid]{gpar}}\code{(lty=2, col = "red")}
#' @param ... Passsed on the other argument in
#' \code{\link[forestplot]{forestplot}}.
#'
#' @return generates a forest plot
#'
#' @seealso \code{\link[forestplot]{forestplot}},
#' \code{\link[forestplot]{fpColors}},
#' \code{\link[forestplot]{fpShapesGp}}, \code{\link[forestplot]{fpLegend}}
#'
#' @examples
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' sce_sub <- allelicRatio(sce_sub)
#' makeForest(sce_sub, showtext = TRUE)
#'
#' # if want to change some properties, like ticks position
#' library(forestplot)
#' xticks <- seq(from = 0, to = 1, by = 0.25)
#' xtlab <- rep(c(TRUE, FALSE), length.out = length(xticks))
#' attr(xticks, "labels") <- xtlab
#' genepoi <- paste0("gene", seq_len(5))
#' ctpoi <- c(1, 3)
#' makeForest(sce_sub, genepoi, ctpoi,
#'   xticks = xticks,
#'   col = fpColors(box = c("blue", "red", "black", "darkgreen"))
#' )
#' @import grid
#' @import forestplot
#' @importFrom utils head
#'
#' @export
makeForest <- function(sce, genepoi, ctpoi = seq_len(nlevels(sce$x)), showtext = FALSE, xticks, boxsize = .25,
                       xlab = "Allelic Ratio", col,
                       grid = structure(seq(0.1, 0.9, 0.1),
                         gp = gpar(lty = 2, col = "#CCCCFF")
                       ), ...) {
  if (missing(xticks)) {
    xticks <- seq(from = 0, to = 1, by = 0.1)
    xtlab <- rep(c(TRUE, FALSE, TRUE, FALSE, FALSE), length.out = length(xticks))
    attr(xticks, "labels") <- xtlab
  }
  if (missing(col)) {
    col <- fpColors(box = brewer.pal(9, "Set1")[seq_len(length(ctpoi))])
  }
  ar <- rowData(sce)[, c(grep("svalue", colnames(rowData(sce)), value = TRUE))] %>% `colnames<-`(levels(sce$x))
  if (showtext) {
    forest_text <- data.frame(`Gene` = rownames(sce), ar[, ctpoi])
  } else {
    forest_text <- data.frame(Gene = rownames(sce)) %>% `rownames<-`(rownames(ar))
  }
  if (missing(genepoi)) {
    smin <- apply(apply(ar, 2, as.numeric), 1, min)
    genepoi <- head(order(smin), min(40, nrow(sce)))
  }
  forest_text <- rbind(colnames(forest_text), forest_text[genepoi, ] %>% as.data.frame())
  message("svalue shown in columns per cell type")
  forestplot::forestplot(forest_text,
    is.summary = c(TRUE, rep(FALSE, nrow(forest_text) - 1)),
    hrzl_lines = list("2" = gpar(lty = 2)),
    line.margin = .1, # We need to add this to avoid crowding
    legend = levels(sce$x)[ctpoi], grid = grid,
    txt_gp = fpTxtGp(
      label = list(gpar(cex = 0.8), gpar(cex = 0.8, col = "#703C3C")),
      ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.9)
    ),
    boxsize = boxsize, graphwidth = unit(15, "cm"),
    mean = rbind(rep(NA, length(ctpoi)), rowData(sce)[genepoi, c(grep("ar", colnames(rowData(sce)), value = TRUE))[ctpoi]] %>% as.matrix()),
    lower = rbind(rep(NA, length(ctpoi)), rowData(sce)[genepoi, c(grep("lower", colnames(rowData(sce)), value = TRUE))[ctpoi]] %>% as.matrix()),
    upper = rbind(rep(NA, length(ctpoi)), rowData(sce)[genepoi, c(grep("upper", colnames(rowData(sce)), value = TRUE))[ctpoi]] %>% as.matrix()),
    clip = c(0.01, 1), xticks = xticks, ref = 0.5,
    col = col, xlab = xlab
  )
}
