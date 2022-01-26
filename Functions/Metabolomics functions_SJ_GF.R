library(EnhancedVolcano)
library(dplyr)
library(forecast)
library(UniprotR)
# library(pathfindR)

BoxCoxTransf = function(x){
  lambda = BoxCox.lambda(x)
  if(lambda == 1){y = log(x)} else {
    y = x^(lambda-1)/lambda}
  if(cor(x,y, use = "pairwise.complete.obs") >= 0){return(y)} else {return(-y)}}

MWAS = function(x,y=group,set1,set2){
  x_new = BoxCoxTransf(x)
  x_new_set1 = x_new[which(y == set1)]
  x_new_set2 = x_new[which(y == set2)]
  p.ttest = t.test(x_new_set1,x_new_set2)$p.val
  fc = mean(x[which(y == set1)]) / mean(x[which(y == set2)])
  # fc = mean(x[which(y == set1)]/x[which(y == set2)])
  log2fc = log(fc,2)
  output = c(log2fc,p.ttest)
  return(output)}

create_results_table = function(df,group1,group2){
  res = t(apply(df,1,MWAS,set1=group1,set2=group2))
  res = data.frame(res)
  colnames(res) = c("log2FC","pval")
  res$FDR = p.adjust(res$pval, method = "fdr") #Default is Holm's
  # Eg. (a<-c(rep(3,4),rep(5,5), rep(1, 11))/100); p.adjust(a); p.adjust(a, "fdr")
  res = res[order(res$pval),]
  return(res)}