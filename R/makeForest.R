#' Plot allelic ratio result as forest
#'
#' Draw a forest plot to visualized cell type specific allelic ratio estimator and confidence interval.
#' It is based on the \pkg{rmeta}-package`s \code{forestplot} function.
#'
#' @param sce A SingleCellExperiment containing colData allelic ratio estimator in the third column
#' and last two column is the confidence interval.
#' @param xticks argument as described in \code{\link[forestplot]{forestplot}}
#' @param boxsize Override the default box size based on precision
#' @param xlab x-axis label. Default is "Allelic Ratio"
#' @param col Set the colors for all the elements. See \code{\link[forestplot]{fpColors}} for details
#' @param grid If you want a discrete gray dashed grid at the level of the ticks
#' you can set this parameter to TRUE.
#' If you set the parameter to a vector of values lines will be drawn at the corresponding positions.
#' If you want to specify the \code{\link[grid]{gpar}} of the lines then either directly pass a \code{\link[grid]{gpar}} object
#' or set the gp attribute e.g. attr(line_vector, "gp") <- \code{\link[grid]{gpar}}(lty=2, col = "red")
#' @param ... Passsed on the other argument in \code{\link[forestplot]{forestplot}}.
#'
#'
#' @seealso \code{\link[forestplot]{forestplot}}, \code{\link[forestplot]{fpColors}},
#' \code{\link[forestplot]{fpShapesGp}}, \code{\link[forestplot]{fpLegend}}
#'
#' @examples
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = 1:4)
#' sce_sub <- wilcoxExt(sce, genecluster = 1)
#' sce_sub <- allelicRatio(sce_sub)
#' makeForest(sce_sub)
#' \dontrun{
#' library(forestplot)
#' xticks <- seq(from = 0, to = 1, by = 0.05)
#' xtlab <- rep(c(TRUE, FALSE), length.out = length(xticks))
#' attr(xticks, "labels") <- xtlab
#'
#' makeForest(sce_sub, xticks, col = fpColors(box = "royalblue"))
#' }
#'
#' @import grid
#' @import forestplot
#'
#' @export
makeForest <- function(sce, xticks, boxsize = 0.1,
                       xlab = "Allelic Ratio", col = fpColors(), grid = structure(c(0.1, 0.5, 0.9),
                         gp = gpar(lty = 2, col = "#CCCCFF")
                       ), ...) {
  if (missing(xticks)) {
    xticks <- seq(from = 0, to = 1, by = 0.05)
    xtlab <- rep(c(TRUE, FALSE), length.out = length(xticks))
    attr(xticks, "labels") <- xtlab
  }

  ar <- apply(metadata(sce)$estimator, 2, as.character)
  forest_text <- rbind(colnames(ar), ar)
  forest_plot <- data.frame(mean = c(NA, ar[, 3]), lower = c(NA, ar[, 4]), upper = c(NA, ar[, 5]))

  forestplot::forestplot(forest_text,
    forest_plot,
    new_page = TRUE, boxsize = boxsize,
    hrzl_lines = list("2" = gpar(lty = 2)), lwd.ci = 2,
    clip = c(0, 1.1), xticks = xticks, grid = grid,
    col = col, xlab = xlab, ...
  )
}
