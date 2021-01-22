fusedlasso<-function(formula,pen.weights,model="quasibinomial",data,lambda="cv1se.dev",k=5,adj.matrix,lambda.length=50L,...){
  misspoi<-which(!is.na(data$ratio))
  # Default is empty list
  if (missing(adj.matrix)) {
    adj.matrix <- list()
  }
  if(model=="quasibinomial"){
    # need to use tryCatch to avoid lambda.max errors
    try1 <- tryCatch({
      fit <- glmsmurf(formula=formula, family=quasibinomial(link = "logit"), data=data,adj.matrix=adj.matrix,
                      weights=data$cts[misspoi], pen.weights="glm.stand", lambda=lambda, 
                      control=list(lambda.length=lambda.length, k=k,...));
      TRUE
    }, error=function(e) {
      message("Failed determining the maximum of lambda, run gaussian model instead")
      fit <-glmsmurf(formula=formula, family=gaussian(), data=data,adj.matrix=adj.matrix,
                     pen.weights="gam.stand", lambda=lambda,control=list(k=k,...))});
  }
  if(model=="gaussian"){
    # need to use tryCatch to avoid lambda.max errors
    try1 <- tryCatch({
       fit <- glmsmurf(formula=formula, family=gaussian(), data=data,
                        pen.weights="glm.stand", lambda=lambda, adj.matrix=adj.matrix,
                        control=list(k=k,lambda.length=lambda.length,...));
       TRUE
    }, error=function(e) {
      message("Failed determining the maximum of lambda, using standardized adaptive GAM weight instead")
      fit <-glmsmurf(formula=formula, family=gaussian(), data=data,adj.matrix=adj.matrix,
                     pen.weights="gam.stand", lambda=lambda,control=list(k=k,lambda.length=lambda.length,...))});
  }
  return(fit)
}

boot_fusedlasso<-function(formula,data,indices,model, lambda1,lambda2,...){
  data_b<-data[indices,]
  poi<-which(!is.na(data_b$ratio))
  if(model=="binomial"){
      fit1 <- glmsmurf(formula=formula, family=binomial(link="logit"), data=data_b,
                      weights=data_b$cts[poi], pen.weights="glm.stand", lambda=lambda1,...);
      co <- coef_reest(fit1)
      co <- co + c(0,rep(co[1],nct-1))
      co <-1/(1+exp(-co))
  }
  if(model=="gaussian"){
    fit2 <- glmsmurf(formula=formula, family=gaussian(), data=data_b,
                    pen.weights="glm.stand", lambda=lambda2,...)
    co <- coef_reest(fit2)
    co <- co + c(0,rep(co[1],nct-1))
  }
  return(co)
}

## Wilcoxon Rank Sum test ###################################
wilcox<-function(data,nct,threshold=0.05,p.adjust.method="none",...){
res <- pairwise.wilcox.test(data$ratio,data$x,p.adjust.method=p.adjust.method,...)
# res <- pairwise.wilcox.test(dat$ratio,dat$x, p.adjust.method ="holm")
adj<-as.data.frame(res$p.value)[lower.tri(res$p.value, diag = T)]
adj<-ifelse(is.nan(adj),1,adj) # Wilcoxon output Nan if ratio of two cell types are exactly same
b<- matrix(0, nct, nct)
b[lower.tri(b, diag=FALSE)]=adj
b2<-b+t(b)
diag(b2)<-1
bb<-ifelse(b2<threshold,1,0) #Binarize p-value to be seen as dismilarity matrix
clust<-hclust(as.dist(bb))
# plot(clust)
my.clusters<-cutree(clust,h=0)
return(my.clusters)
}


library(tidyverse)
wilcox_adj<-function(data,nct,k,threshold,p.adjust.method="none",...){
  set.seed(as.numeric(Sys.Date()))
  # foldInds <- createFolds(data$x, k=k, list=TRUE, returnTrain=T)
  out<-list()
  obj<-sapply (1:length(threshold), function(j){
    # sae<-sapply(1:k,function(i){
    #   fit3<-wilcox(data[foldInds[[i]],],nct=nct,threshold=threshold[j])
    #   label<-tibble(type=factor(seq_along(1:nct)),par=fit3)
    #   test<-data[-foldInds[[i]],]
    #   test2<-test %>% left_join(label,by=c("x"="type"))
    #   test2<-test2 %>% group_by(par) %>% mutate(grpmean=mean(ratio,na.rm = T))
    #   a<-sum((test2$ratio-test2$grpmean)^2,na.rm = T) #se
    #   b<-sum(abs(test2$ratio-test2$grpmean),na.rm = T) #ae
    #   return(c(a,b))
    # })
    fit<-wilcox(data,nct=nct,p.adjust.method=p.adjust.method,threshold=threshold[j],...)
    label<-tibble(type=factor(seq_along(1:nct)),par=fit)
    data2<-data %>% left_join(label,by=c("x"="type"))
    data2<-data2 %>% group_by(par) %>% mutate(grpmean=mean(ratio,na.rm = T))
    loss1<-nrow(data)*log(sum((data2$ratio-data2$grpmean)^2,na.rm = T)/nrow(data2))+length(unique(fit))*log(nrow(data2)) #BIC

    ## CV1mse
    # mean<-mean(sae[1,])
    # sd<-sd(sae[1,])/sqrt(k)
    # absolute.mean<-mean(sae[2,])
    # absolute.sd<-sd(sae[2,])/sqrt(k)
    out[["cl"]]<-fit
    out[["loss1"]]<-loss1
    # out[["smean"]]<-mean
    # out[["ssd"]]<-sd
    # out[["amean"]]<-absolute.mean
    # out[["asd"]]<-absolute.sd
    return(out)
  }) 
  cl<-as_tibble(do.call(rbind, obj[seq(1,length(obj), by = 2)]))
  loss1<-as_tibble(do.call(rbind, obj[seq(2,length(obj), by = 2)]))
  # smean<-as_tibble(do.call(rbind, obj[seq(3,length(obj), by = 6)]))
  # ssd<-as_tibble(do.call(rbind, obj[seq(4,length(obj), by = 6)]))
  # amean<-as_tibble(do.call(rbind, obj[seq(5,length(obj), by = 6)]))
  # asd<-as_tibble(do.call(rbind, obj[seq(6,length(obj), by = 6)]))
  return(cl[which.min(unlist(loss1)),])
  # return(cl<-list(cl[which.min(unlist(loss1)),],
  #         cl[which(smean < min(smean) + mean(as.numeric(unlist(ssd)))/sqrt(k))[1],],
  #        cl[which(amean < min(amean) + mean(as.numeric(unlist(asd)))/sqrt(k))[1],]))
}


