context("fusedLasso")
library(S4Vectors)

test_that("check inputs", {

  # missing partition step
  sce <- makeSimulatedData()
  sce <- preprocess(sce)
  sce <- geneCluster(sce, G = seq_len(4))
  expect_error(allelicRatio(sce))

  # missing input R
  sce_sub <- wilcoxExt(sce,genecluster = 1)
  expect_error(allelicRatio(sce_sub,method="bootstrap"), "bootstrap replicates")

})

test_that("basic betabinom analyses", {

  sce <- makeSimulatedData()
  sce <- preprocess(sce)
  sce <- geneCluster(sce, G = seq_len(4))
  sce_sub <- wilcoxExt(sce,genecluster = 1)
  sce_sub <- allelicRatio(sce_sub)
  expect_true(!is.null(metadata(sce_sub)$estimator))
})

test_that("bootstrap CI", {

  sce <- makeSimulatedData()
  sce <- preprocess(sce)
  sce <- geneCluster(sce, G = seq_len(4))
  sce_sub <- wilcoxExt(sce,genecluster = 1)
  sce_sub <- allelicRatio(sce_sub,model="bootstrap", R=5)
  expect_true(!is.null(metadata(sce_sub)$estimator))
})
