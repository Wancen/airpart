#' @importFrom smurf glmsmurf
#' @export
boot_fusedlasso<-function(formula,data,indices,model, lambda,...){
  data_b<-data[indices,]
  poi<-which(!is.na(data_b$ratio))
  if(model=="binomial"){
    fit1 <- smurf::glmsmurf(formula=formula, family=binomial(link="logit"), data=data_b,
                     weights=data_b$cts[poi], pen.weights="glm.stand", lambda=lambda1,...);
    co <- coef_reest(fit1)
    co <- co + c(0,rep(co[1],nct-1))
    co <-1/(1+exp(-co))
  }
  if(model=="gaussian"){
    fit2 <- smurf::glmsmurf(formula=formula, family=gaussian(), data=data_b,
                     pen.weights="glm.stand", lambda=lambda2,...)
    co <- coef_reest(fit2)
    co <- co + c(0,rep(co[1],nct-1))
  }
  return(co)
}
