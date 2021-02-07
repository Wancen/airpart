#' @import smurf
#'
#' @export
#' @title Run generalized fused lasso to part cell-type-specific allelic imbalance across all cell types
#'
#' @description Fit generalized fused lasso with either binomial(link="logit") or gaussian likelihood.
#'
#' @param formula A \code{\link[stats]{formula}} object describing the model to be fitted.
#' Penalties are specified using the \code{\link{p}} function.
#' @param model Either "binomial" or "gaussian" used to model the generalized fused lasso
#' @param data A data frame containing the model response and predictors for \code{N} observations
#' @param pen.weights Either a string describing the method to compute the penalty weights:
#'  \itemize{
#'                  \item \code{"eq"} (equal penalty weights),
#'                  \item \code{"stand"} (standardization penalty weights),
#'                  \item \code{"glm"} (adaptive GLM penalty weights),
#'                  \item \code{"glm.stand"} (default; stand. ad. GLM penalty weights),
#'                  \item \code{"gam"} (ad. GAM penalty weights),
#'                  \item \code{"gam.stand"} (stand. ad. GAM penalty weights);
#'               }
#'                or a list with the penalty weight vector per predictor. This list should have length equal to the number of predictors and predictor names as element names.
#' @param lambda Either the penalty parameter, a positive number; or a string describing the method and measure used to select the penalty parameter:
#'               \itemize{
#'                  \item \code{"is.aic"} (in-sample; Akaike Information Criterion (AIC)),
#'                  \item \code{"is.bic"} (in-sample; Bayesian Information Criterion (BIC)),
#'                  \item \code{"is.gcv"} (in-sample; Generalized Cross-Validation (GCV) score),
#'                  \item \code{"oos.dev"} (out-of-sample; deviance),
#'                  \item \code{"oos.mse"} (out-of-sample; Mean Squared Error (MSE)),
#'                  \item \code{"oos.dss"} (out-of-sample; Dawid-Sebastiani Score (DSS)),
#'                  \item \code{"cv.dev"} (cross-validation (CV); deviance),
#'                  \item \code{"cv.mse"} (CV; MSE),
#'                  \item \code{"cv.dss"} (CV; DSS),
#'                  \item \code{"cv1se.dev"} (default; CV with one standard error (SE) rule; deviance),
#'                  \item \code{"cv1se.mse"} (CV with one SE rule; MSE),
#'                  \item \code{"cv1se.dss"} (CV with one SE rule; DSS).
#'               }
#'               E.g. \code{"is.aic"} indicates in-sample selection of lambda with the AIC as measure.
#'               When \code{lambda} is missing or \code{NULL}, it will be selected using cross-validation with the one standard error rule and the deviance as measure (\code{"cv1se.dev"}).
#' @param niter number of iteration to run, recommend run 5 times if allelic ratio difference within [0.05,0.1]
#' @param adj.matrix A named list containing the adjacency matrices (a.k.a. neighbor matrices) for each of the predictors with a Graph-Guided Fused Lasso penalty.
#'                The list elements should have the names of the corresponding predictors. If only one predictor has a Graph-Guided Fused Lasso penalty,
#'                it is also possible to only give the adjacency matrix itself (not in a list).
#'
#' @return An object of class 'glmsmurf' is returned.
#'
#' @details See the package vignette for more details and a complete description of a use case.
#'
#' @seealso \code{\link[smurf]{glmsmurf}}, \code{\link[smurf]{glmsmurf.control}}, \code{\link[smurf]{p}}, \code{\link[stats]{glm}}
fusedlasso<-function(formula,model="binomial",data,pen.weights,lambda="cv1se.dev",k=5,niter=1,adj.matrix,lambda.length=25L,...){
  misspoi<-which(!is.na(data$ratio))
  nct<-length(levels(data$x))
  # Default is empty list
  if (missing(adj.matrix)) {
    adj.matrix <- list()
  }
  if(model=="binomial"){
    # need to use tryCatch to avoid lambda.max errors
    try1 <- tryCatch({

      coef<-sapply(1:niter,function(t){
        fit<-smurf::glmsmurf(formula=formula, family=binomial(link = "logit"), data=data,adj.matrix=adj.matrix,
                             weights=data$cts[misspoi], pen.weights="glm.stand", lambda=lambda,
                             control=list(lambda.length=lambda.length, k=k,...));
        co <- coef_reest(fit)
        co <- co + c(0,rep(co[1],nct-1))
        if(nct<10){
          index.1<-which(rowMeans(fit$lambda.measures$dev) < min(rowMeans(fit$lambda.measures$dev)) +0.5*mean(matrixStats::rowSds(fit$lambda.measures$dev)/sqrt(5)))[1]
          fit2<-smurf::glmsmurf(formula=formula, family=binomial(link = "logit"), data=data,adj.matrix=adj.matrix,
                                weights=data$cts[misspoi], pen.weights="glm.stand",lambda=fit$lambda.vector[index.1],
                                control = list(...)) #0.5 SE
          co <- coef_reest(fit2)
          co <- co + c(0,rep(co[1],nct-1))
        }
        return(co)
      })
      TRUE
    }, error=function(e) {
      message("Failed determining the maximum of lambda, try run other weights or gaussian model instead")});
  }
  if(model=="gaussian"){
    # need to use tryCatch to avoid lambda.max errors
    try1 <- tryCatch({
      coef<-sapply(1:niter,function(t){
        fit<-smurf::glmsmurf(formula=formula, family=gaussian(), data=data,adj.matrix=adj.matrix,
                             weights=data$cts[misspoi], pen.weights="glm.stand", lambda=lambda,
                             control=list(lambda.length=lambda.length, k=k,...));
        co <- coef_reest(fit)
        co <- co + c(0,rep(co[1],nct-1))
        if(nct<10){
          index.1<-which(rowMeans(fit$lambda.measures$dev) < min(rowMeans(fit$lambda.measures$dev)) +0.5*mean(matrixStats::rowSds(fit$lambda.measures$dev)/sqrt(5)))[1]
          fit2<-smurf::glmsmurf(formula=formula, family=gaussian(), data=data,adj.matrix=adj.matrix,
                                weights=data$cts[misspoi], pen.weights="glm.stand",lambda=fit$lambda.vector[index.1]) #0.5 SE
          co <- coef_reest(fit2)
          co <- co + c(0,rep(co[1],nct-1))
        }
        return(co)
      })
       TRUE
    }, error=function(e) {
      message("Failed determining the maximum of lambda, try run other weights instead")});

  }
  cl <-apply(coef, 2, function(x) fmatch(x,unique(x)))
  return(cl)
}





