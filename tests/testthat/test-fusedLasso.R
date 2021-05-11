library(S4Vectors)

test_that("wrong inputs", {

  # missed input arrays
  set.seed(2021)
  sce <- makeSimulatedData()
  expect_error(fusedLasso(sce, ncores = 2))

  # no genecluster identified
  sce <- preprocess(sce)
  expect_error(fusedLasso(sce, ncores = 2), "No gene cluster number")

  # no cluster column
  expect_error(fusedLasso(sce, model = "binomial", genecluster = 1, ncores = 2))

  # wrong cell state column names
  names(colData(sce)) <- "ct"
  expect_error(fusedLasso(sce, model = "binomial", genecluster = 1, ncores = 2))
})

test_that("basic fusedLasso analyses", {
  set.seed(2021)
  sce <- makeSimulatedData()
  sce <- preprocess(sce)
  sce <- geneCluster(sce, G = seq_len(4), plot = FALSE)
  sce_sub <- fusedLasso(sce, model = "binomial", genecluster = 1, ncores = 2)
  expect_true(is.numeric(metadata(sce_sub)$lambda))
})
