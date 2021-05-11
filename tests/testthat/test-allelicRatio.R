library(S4Vectors)

test_that("check inputs", {

  # missing partition step
  set.seed(2021)
  sce <- makeSimulatedData()
  sce <- preprocess(sce)
  sce <- geneCluster(sce, G = seq_len(4))
  expect_error(allelicRatio(sce))
})

test_that("basic betabinom analyses", {
  set.seed(2021)
  sce <- makeSimulatedData()
  sce <- preprocess(sce)
  sce <- geneCluster(sce, G = seq_len(4))
  sce_sub <- wilcoxExt(sce, genecluster = 1)
  sce_sub <- allelicRatio(sce_sub)
  expect_true(!is.null(rowData(sce_sub)[, c(grep("ar", colnames(rowData(sce_sub)), value = TRUE))]))
})
