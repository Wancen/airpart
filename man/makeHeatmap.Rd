% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeHeatmap.R
\name{makeHeatmap}
\alias{makeHeatmap}
\title{Plot allelic ratio as heatmap}
\usage{
makeHeatmap(
  sce,
  assay = c("ratio_pseudo", "ratio", "counts"),
  genecluster = NULL,
  show_row_names = FALSE,
  order_by_group = TRUE,
  ...
)
}
\arguments{
\item{sce}{SingleCellExperiment}

\item{assay}{the assay to be plotted. Choices are \code{"ratio_pseudo"} which is the default,
\code{"ratio"}, \code{"counts"}.}

\item{genecluster}{an integer indicates which gene cluster heatmap want to be returned.}

\item{show_row_names}{show row names or not}

\item{order_by_group}{indicate whether order by group or order by cell types}

\item{...}{Passsed on the other argument in
\code{\link[ComplexHeatmap]{Heatmap}}.}
}
\value{
generates a heatmap
}
\description{
Plot allelic ratio as heatmap
}
\examples{
set.seed(2021)
sce <- makeSimulatedData(p.vec = c(0.3, 0.5, 0.5, 0.3), ncl = 1)
sce <- preprocess(sce)
# display allelic ratio pattern in whole dataset
makeHeatmap(sce)

sce <- geneCluster(sce, G = seq_len(4), plot = FALSE)
sce_sub <- wilcoxExt(sce, genecluster = 1)
# display specific gene cluster partition result
makeHeatmap(sce_sub)
# display by cell type orders
makeHeatmap(sce_sub, order_by_group = FALSE)
}
