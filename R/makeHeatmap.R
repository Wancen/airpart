#' Plot allelic ratio as heatmap
#'
#' @param se SummarizedExpeirment
#'
#' @importFrom ComplexHeatmap Heatmap
#' @export
makeRatioHeatmap <- function(se,show_row_names = F) {
  m<-assays(se)[["ratio_pseudo"]]
  if(is.null(se$part)){
    ha=HeatmapAnnotation(`cell type`=se$x,border = F,col = list(`cell type` = structure(brewer.pal(nlevels(se$x), "Set3"),names = levels(se$x))))
    Heatmap(m, name = "Allelic Ratio", cluster_columns = F,cluster_rows = F,show_column_names = F,show_row_names = show_row_names,top_annotation=ha)
  }else{
    split = se$part
    ha = HeatmapAnnotation(group = anno_block(gp = gpar(fill = brewer.pal(9, "Pastel1")[1:nlevels(split)]),
                                              labels_gp = gpar(col = "white", fontsize = 10),
                                              labels = paste0("group",1:nlevels(split))),
                           `cell type`=se$x,border = F,
                           col = list(`cell type` = structure(brewer.pal(nlevels(se$x), "Set3"),names = levels(se$x))))
    Heatmap(m, name = "Allelic Ratio", column_split = split,cluster_columns = F,cluster_rows = F,show_column_names = F,show_row_names = show_row_names,top_annotation=ha)
  }
}





