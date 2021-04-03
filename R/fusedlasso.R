#' Generalized fused lasso to partition cell types by allelic imbalance
#'
#' Fits generalized fused lasso with either binomial(link="logit")
#' or gaussian likelihood, leveraging functions from the
#' \code{smurf} package.
#'
#' @param sce A SingleCellExperiment containing assays (\code{"ratio"},
#' \code{"counts"}) and colData \code{"x"}
#' @param formula A \code{\link[stats]{formula}} object which will typically be
#' fused lasso penalty:
#' \code{ratio ~ p(x, pen="gflasso")}. Another possibility would be to use
#' the Graph-Guided Fused Lasso penalty:
#' \code{f <- ratio ~ p(x, pen = "ggflasso")}
#' See \code{\link[smurf]{glmsmurf}} for more details
#' @param model Either \code{"binomial"} or \code{"gaussian"} used to fit
#' the generalized fused lasso
#' @param genecluster which gene cluster result want to be returned.
#' Usually identified interesting gene cluster pattern by
#' \code{\link{summaryAllelicRatio}}
#' @param niter number of iteration to run; recommended to run 5 times
#' if allelic ratio differences are within [0.05,0.1]
#' @param pen.weights argument as described in \code{\link[smurf]{glmsmurf}}
#' @param lambda argument as described in \code{\link[smurf]{glmsmurf}}.
#' Default lambda is determined by \code{"cv1se.dev"}
#' (cross validation within 1 standard error rule(SE); deviance)
#' @param k number of cross-validation folds
#' @param adj.matrix argument as described in \code{\link[smurf]{glmsmurf}}
#' @param lambda.length argument as described in \code{\link[smurf]{glmsmurf}}
#' @param se.rule.nct the number of cell types to trigger
#' another SE based rule
#' (to prioritize larger models, less fusing,
#' good to detect 0.05 allelic ratio difference).
#' When the number of cell types is less than or equal to this value, the
#' \code{se.rule.mult} SE rule is used
#' @param se.rule.mult the multiplier of the SE in determining the lambda:
#' the chosen lambda is within \code{se.rule.mult} x SE of the minimum
#' deviance. Default is 0.5 SE
#' @param ... additional arguments passed to \code{\link[smurf]{glmsmurf}}
#'
#' @return A matrix grouping factor partition
#' and the penalized parameter lambda
#' are returned in metadata \code{"partition"} and \code{"lambda"}.
#' Partation also stored in colData\code{"part"}.
#'
#' @details Usually, we used a Generalized Fused Lasso penalty for the
#' cell states in order to regularize all possible coefficient differences.
#' Another possibility would be to use the Graph-Guided Fused Lasso penalty
#' to only regularize the differences of coefficients of neighboring
#' cell states.
#' When using a Graph-Guided Fused Lasso penalty, the adjacency matrix
#' corresponding to the graph needs to be provided. The elements of this
#' matrix are zero when two levels are not connected, and one when they are
#' adjacent.
#'
#' See the package vignette for more details and a complete description of a
#' use case.
#'
#' @references
#'
#' This function leverages the glmsmurf function from the smurf package.
#' For more details see the following manuscript:
#'
#' Devriendt S, Antonio K, Reynkens T, et al.
#' Sparse regression with multi-type regularized feature modeling[J].
#' Insurance: Mathematics and Economics, 2021, 96: 248-261.
#'
#' @seealso \code{\link[smurf]{glmsmurf}},
#' \code{\link[smurf]{glmsmurf.control}},
#' \code{\link[smurf]{p}}, \code{\link[stats]{glm}}
#'
#' @examples
#' library(S4Vectors)
#' library(smurf)
#' sce <- makeSimulatedData()
#' sce <- preprocess(sce)
#' sce <- geneCluster(sce, G = seq_len(4))
#' f <- ratio ~ p(x, pen = "gflasso") # formula for the GFL
#' sce_sub <- fusedLasso(sce,
#'   formula = f, model = "binomial", genecluster = 1,
#'   ncores = 2, se.rule.nct = 3
#' )
#' metadata(sce_sub)$partition
#' metadata(sce_sub)$lambda
#'
#' # Suppose we have 4 cell states, if we don't want cell state 1
#' # to be grouped together with other cell states
#' adj.matrix <- 1 - diag(4)
#' colnames(adj.matrix) <- rownames(adj.matrix) <- levels(sce$x)
#' adj.matrix[1, c(2, 3, 4)] <- 0
#' adj.matrix[c(2, 3, 4), 1] <- 0
#' f <- ratio ~ p(x, pen = "ggflasso") # use graph-guided fused lasso
#' sce_sub <- fusedLasso(sce,
#'   formula = f, model = "binomial", genecluster = 1,
#'   lambda = 0.5, ncores = 2, se.rule.nct = 3,
#'   adj.matrix = adj.matrix
#' )
#' metadata(sce_sub)$partition
#' @import smurf
#' @importFrom matrixStats rowSds
#' @importFrom stats binomial gaussian
#'
#' @export
fusedLasso <- function(sce, formula, model = "binomial",
                       genecluster, niter = 1,
                       pen.weights, lambda = "cv1se.dev", k = 5,
                       adj.matrix, lambda.length = 25L,
                       se.rule.nct = 8,
                       se.rule.mult = 0.5,
                       ...) {
  if (missing(genecluster)) {
    stop("No gene cluster number")
  }

  stopifnot(c("ratio", "counts") %in% assayNames(sce))
  stopifnot("x" %in% names(colData(sce)))
  stopifnot("cluster" %in% names(rowData(sce)))

  if (missing(formula)) {
    formula <- ratio ~ p(x, pen="gflasso")
  }
  # default is empty list
  if (missing(adj.matrix)) {
    adj.matrix <- list()
  }
  sce_sub <- sce[rowData(sce)$cluster == genecluster, ]
  cl_ratio <- as.vector(unlist(assays(sce_sub)[["ratio"]]))
  cl_total <- as.vector(unlist(counts(sce_sub)))
  dat <- data.frame(
    ratio = cl_ratio,
    x = factor(rep(sce_sub$x, each = length(sce_sub))),
    cts = cl_total
  )
  dat <- dat[!is.nan(dat$ratio), ]
  if (model == "binomial") {
    fam <- binomial(link = "logit")
    msg <- "Failed determining max lambda, try other weights or gaussian model"
    weight <- dat$cts
  } else {
    fam <- gaussian()
    msg <- "Failed determining max of lambda, try other weights"
    weight <- NULL
  }
  nct <- nlevels(sce$x)
  # need to use tryCatch to avoid lambda.max errors
  res <- tryCatch(
    {
        vapply(seq_len(niter), function(t) {
          fit <- smurf::glmsmurf(
          formula = formula, family = fam,
          data = dat, adj.matrix = adj.matrix,
          weights = weight,
          pen.weights = "glm.stand", lambda = lambda,
          control = list(lambda.length = lambda.length, k = k, ...)
        )
        co <- coef_reest(fit)
        co <- co + c(0, rep(co[1], nct - 1))
        lambda <- fit$lambda
        # if number of cell types is 'se.rule.nct' or less:
        if (nct <= se.rule.nct) {
          # choose lambda by the lowest deviance within 'se.rule.mult'
          # standard error of the min
          mean.dev <- rowMeans(fit$lambda.measures$dev)
          min.dev <- min(mean.dev)
          sd.dev <- matrixStats::rowSds(fit$lambda.measures$dev)
          se.dev <- mean(sd.dev) / sqrt(k)
          idx <- which(mean.dev < min.dev + se.rule.mult * se.dev)[1]
          # this is faster, running the GFL for a single lambda value
          fit2 <- smurf::glmsmurf(
            formula = formula, family = fam,
            data = dat, adj.matrix = adj.matrix,
            weights = weight, pen.weights = "glm.stand",
            lambda = fit$lambda.vector[idx],
            control = list(...)
          )
          # rearrange coefficients so not comparing to reference cell type
          co <- coef_reest(fit2)
          co <- co + c(0, rep(co[1], nct - 1))
          lambda <- fit2$lambda
        }
        return(c(co, lambda))
      }, double(nct+1))
    },
    error = function(e) {
      message(msg)
      return(NA)
    }
  )
  if (length(res) == 1 && is.na(res))
    stop("Error occurred in attempting to run fused lasso")
  if (niter == 1) {
    coef <- res[seq_len(nct), ]
    lambda <- unname(res[nct + 1, ])
    part <- match(coef, unique(coef)) %>% as.factor()
  } else {
    # multiple partitions
    coef <- res[seq_len(nct), ]
    lambda <- res[nct + 1, ]
    part <- apply(coef, 2, function(z) match(z, unique(z)))
    colnames(part) <- paste0("part", seq_len(niter))
    names(lambda) <- paste0("part", seq_len(niter))
  }
  cl <- data.frame(part, x = levels(sce_sub$x))
  cd <- colData(sce_sub)
  cd2 <- cd[, !names(cd) %in% c("part", "rowname")] %>%
    as.data.frame() %>%
    setNames(names(cd)[!names(cd) %in% c("part", "rowname")])
  coldata <- DataFrame(rowname = colnames(sce_sub), cd2)
  coldata <- merge(coldata, cl, by = "x", sort = FALSE) %>%
    DataFrame()
  rownames(coldata) <- coldata$rowname
  colData(sce_sub) <- coldata
  metadata(sce_sub)$partition <- cl
  metadata(sce_sub)$lambda <- lambda
  return(sce_sub)
}
