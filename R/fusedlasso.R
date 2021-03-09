#' Generalized fused lasso to partition cell types by allelic imbalance
#'
#' Fits generalized fused lasso with either binomial(link="logit")
#' or gaussian likelihood, leveraging functions from the
#' \code{smurf} package.
#'
#' @param formula A \code{\link[stats]{formula}} object which will typically be:
#' \code{ratio ~ p(x, pen="gflasso")}.
#' See \code{\link[smurf]{glmsmurf}} for more details
#' @param model Either \code{"binomial"} or \code{"gaussian"} used to fit the generalized fused lasso
#' @param se A SummarizedExpeirment containing assays (\code{"ratio"},
#' \code{"total"}) and colData \code{"x"}
#' @param genecluster which gene cluster result want to be returned.
#' Usually identified interesting gene cluster pattern by \code{\link{summaryAllelicRatio}}
#' @param niter number of iteration to run; recommended to run 5 times
#' if allelic ratio differences are within [0.05,0.1]
#' @param pen.weights argument as described in \code{\link[smurf]{glmsmurf}}
#' @param lambda argument as described in \code{\link[smurf]{glmsmurf}}
#' @param k number of cross-validation folds
#' @param adj.matrix argument as described in \code{\link[smurf]{glmsmurf}}
#' @param lambda.length argument as described in \code{\link[smurf]{glmsmurf}}
#' @param se.rule.nct the number of cell types to trigger an SE based rule
#' (to prioritize smaller models, more fusing). When the number of cell types
#' is less than or equal to this value, the SE rule is used
#' @param se.rule.mult the multiplier of the SE in determining the lambda:
#' the chosen lambda is within \code{se.rule.mult} x SE of the minimum
#' deviance
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
fusedLasso <- function(formula, model="binomial", se, genecluster, niter=1,
                       pen.weights, lambda="cv1se.dev", k=5,
                       adj.matrix, lambda.length=25L,
                       se.rule.nct=8,
                       se.rule.mult=0.5,
                       ...) {
  if (missing(genecluster)) {
    stop("No gene cluster number")
  }
  stopifnot(c("ratio","total") %in% assayNames(se))

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
  se_sub<-se[rowData(se)$cluster == genecluster, ]
  cl_ratio <- as.vector(unlist(assays(se_sub)[["ratio"]]))
  cl_total <- as.vector(unlist(assays(se_sub)[["total"]]))
  dat <- data.frame(ratio=cl_ratio,
                    x=factor(rep(se_sub$x,each=length(se_sub))),
                    cts=cl_total)
  dat <- dat[!is.nan(dat$ratio),]
  nct<-nlevels(se$x)
  # need to use tryCatch to avoid lambda.max errors
  try <- tryCatch({
    coef <- sapply(1:niter, function(t) {
      fit <- smurf::glmsmurf(formula=formula, family=fam,
                             data=dat, adj.matrix=adj.matrix,
                             weights=dat$cts,
                             pen.weights="glm.stand", lambda=lambda,
                             control=list(lambda.length=lambda.length, k=k, ...))
      co <- coef_reest(fit)
      co <- co + c(0,rep(co[1],nct-1))
      # if number of cell types is 'se.rule.nct' or less:
      if (nct <= se.rule.nct) {
        # choose lambda by the lowest deviance within 'se.rule.mult' standard error of the min
        mean.dev <- rowMeans(fit$lambda.measures$dev)
        min.dev <- min(mean.dev)
        sd.dev <- matrixStats::rowSds(fit$lambda.measures$dev)
        se.dev <- mean(sd.dev)/sqrt(k)
        idx <- which(mean.dev < min.dev + se.rule.mult * se.dev)[1]
        # this is faster, running the GFL for a single lambda value
        fit2 <- smurf::glmsmurf(formula=formula, family=fam,
                                data=dat, adj.matrix=adj.matrix,
                                weights=dat$cts, pen.weights="glm.stand",
                                lambda=fit$lambda.vector[idx],
                                control = list(...))
        # rearrange coefficients so not comparing to reference cell type
        co <- coef_reest(fit2)
        co <- co + c(0,rep(co[1],nct-1))
      }
      return(co)
    })
    TRUE
  }, error=function(e) {
    message(msg)
  })
  if (niter == 1) {
    part <- match(coef[,1], unique(coef[,1]))
  } else {
    # multiple partitions
    part <- apply(coef, 2, function(z) match(z, unique(z)))
    colnames(part) <- paste0("part",seq_len(niter))
  }
  cl <- data.frame(part,x =levels(se_sub$x))
  colData(se_sub)<-merge(colData(se_sub),cl,by="x")%>% DataFrame() %>% `row.names<-`(colnames(se_sub))
  metadata(se_sub)$partition<-cl
  return(se_sub)
}
