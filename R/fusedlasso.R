#' Generalized fused lasso to partition cell types by allelic imbalance
#'
#' Fits generalized fused lasso with either binomial(link="logit")
#' or gaussian likelihood, leveraging functions from the
#' \code{smurf} package.
#'
#' @param formula A \code{\link[stats]{formula}} object describing the model to be fitted.
#' See \code{\link[smurf]{glmsmurf}} for more details
#' @param model Either \code{"binomial"} or \code{"gaussian"} used to fit the generalized fused lasso
#' @param data A data frame containing the allelic ratio, the cell assignment,
#' and the total counts for weighting (total counts should be specified as a
#' column \code{cts}). The allelic ratio and cell assignment are as named in
#' \code{formula}, e.g. \code{ratio ~ p(x, pen="gflasso")}
#' @param pen.weights argument as described in \code{\link[smurf]{glmsmurf}}
#' @param lambda argument as described in \code{\link[smurf]{glmsmurf}}
#' @param niter number of iteration to run; recommended to run 5 times
#' if allelic ratio differences are within [0.05,0.1]
#' @param adj.matrix argument as described in \code{\link[smurf]{glmsmurf}}
#' @param ... additional arguments passed to \code{\link[smurf]{glmsmurf}}
#'
#' @return An object of class 'glmsmurf' is returned.
#'
#' @details See the package vignette for more details and a complete description of a use case.
#'
#' @references
#'
#' This function leverages the glmsmurf function from the smurf package.
#' For more details see the following manuscript:
#'
#' Sander Devriendt, Katrien Antonio, Tom Reynkens, Roel Verbelen
#' "Sparse Regression with Multi-type Regularized Feature Modeling"
#' 2018. arXiv:1810.03136 [stat.CO].
#' 
#' @seealso \code{\link[smurf]{glmsmurf}}, \code{\link[smurf]{glmsmurf.control}},
#' \code{\link[smurf]{p}}, \code{\link[stats]{glm}}
#'
#' @import smurf
#' @importFrom matrixStats rowSds
#'
#' @export
fusedlasso <- function(formula, model="binomial", data,
                       pen.weights, lambda="cv1se.dev",
                       k=5, niter=1, adj.matrix,
                       lambda.length=25L, ...) {
  # TODO: can we just use this?
  data <- data[!is.na(data$ratio),]
  nct <- length(levels(data$x))
  # default is empty list
  if (missing(adj.matrix)) {
    adj.matrix <- list()
  }
  if (model=="binomial") {
    fam <- binomial(link = "logit")
    msg <- "Failed determining max of lambda, try other weights or gaussian model"
  } else {
    fam <- gaussian()
    msg <- "Failed determining max of lambda, try other weights"
  }
  # need to use tryCatch to avoid lambda.max errors
  try <- tryCatch({
    coef <- sapply(1:niter, function(t) {
      fit <- smurf::glmsmurf(formula=formula, family=fam,
                             data=data, adj.matrix=adj.matrix,
                             weights=data$cts,
                             pen.weights="glm.stand", lambda=lambda,
                             control=list(lambda.length=lambda.length, k=k,...))
      co <- coef_reest(fit)
      co <- co + c(0,rep(co[1],nct-1))
      if (nct <= 8) {
        # TODO: can you put comments here what this code is for?
        # probably it should also be described in the details above
        idx <- which(rowMeans(fit$lambda.measures$dev) < min(rowMeans(fit$lambda.measures$dev)) +
                         0.5 * mean(matrixStats::rowSds(fit$lambda.measures$dev)/sqrt(5)))[1]
        fit2 <- smurf::glmsmurf(formula=formula, family=fam,
                                data=data, adj.matrix=adj.matrix,
                                weights=data$cts, pen.weights="glm.stand",
                                lambda=fit$lambda.vector[idx],
                                control = list(...)) # 0.5 SE
        co <- coef_reest(fit2)
        co <- co + c(0,rep(co[1],nct-1))
      }
      return(co)
    })
    TRUE
  }, error=function(e) {
    message(msg)
  })
  cl <- apply(coef, 2, function(x) match(x,unique(x)))
  rownames(cl) <- levels(x)
  return(cl)
}
