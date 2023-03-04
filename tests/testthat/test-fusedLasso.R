library(S4Vectors)

test_that("wrong inputs", {

  # missed input arrays
  set.seed(2021)
  sce <- makeSimulatedData()
  expect_error(fusedLasso(sce, ncores = 1))

  # no genecluster identified
  sce <- preprocess(sce)
  expect_error(fusedLasso(sce, ncores = 1), "No gene cluster number")

  # no cluster column
  expect_error(fusedLasso(sce, model = "binomial", genecluster = 1, ncores = 1))

  # wrong cell state column names
  names(colData(sce)) <- "ct"
  expect_error(fusedLasso(sce, model = "binomial", genecluster = 1, ncores = 1))
})

test_that("basic fusedLasso analyses", {
  set.seed(2021)
  sce <- makeSimulatedData()
  sce <- preprocess(sce)
  sce <- geneCluster(sce, G = seq_len(4), plot = FALSE)
  sce_sub <- fusedLasso(sce, model = "binomial", genecluster = 1, ncores = 1)
  expect_true(is.numeric(metadata(sce_sub)$lambda))

  # gene and donor in formula
  colData(sce)$donor <- factor(rep(1:12, each=10))
  f <- ratio ~ p(x, pen = "gflasso") + gene + donor
  sce_sub <- fusedLasso(sce[1:5,], formula=f,
                        model = "binomial",
                        genecluster = 1, ncores = 1)
  expect_true(is.factor(metadata(sce_sub)$partition$part))
})
