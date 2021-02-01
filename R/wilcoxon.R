#' @export
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
  clust<-hclust(as.dist(bb)) #hierarchical cluster on adjacency matrix
  my.clusters<-cutree(clust,h=0)
  return(my.clusters)
}

wilcox_adj<-function(data,nct,k,threshold,p.adjust.method="none",...){
  set.seed(as.numeric(Sys.Date()))
  out<-list()
  obj<-sapply (1:length(threshold), function(j){
    fit<-wilcox(data,nct=nct,p.adjust.method=p.adjust.method,threshold=threshold[j],...)
    label<-tibble(type=factor(seq_along(1:nct)),par=fit)
    data2<-data %>% left_join(label,by=c("x"="type"))
    data2<-data2 %>% group_by(par) %>% mutate(grpmean=mean(ratio,na.rm = T))
    # loss function
    loss1<-nrow(data)*log(sum((data2$ratio-data2$grpmean)^2,na.rm = T)/nrow(data2))+length(unique(fit))*log(nrow(data2))
    out[["cl"]]<-fit
    out[["loss1"]]<-loss1
    return(out)
  })
  cl<-as_tibble(do.call(rbind, obj[seq(1,length(obj), by = 2)]))
  loss1<-as_tibble(do.call(rbind, obj[seq(2,length(obj), by = 2)]))
  return(cl[which.min(unlist(loss1)),])
}

