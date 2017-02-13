require(sva)
require(data.table)
require(ggplot2)
require(glmnet)
require(survival)
#require(reshape2)
require(rms)
require(caret)
require(broom)
#require(affy)
#require(affydata)
options(warn=-1) 
# 
# combat_expr_data = fread('D:/R_workspace/survival/PPT/combat_expr_data.csv')
# patient_info = read.csv('D:/R_workspace/survival/PPT/patient_info.csv',row.names = 1)
# pvalue_logrank = read.csv('D:/R_workspace/survival/PPT/pvalue_logrank.csv',row.names = 1)
# 
# combat_expr_data = as.data.frame(combat_expr_data)
# rownames(combat_expr_data) = combat_expr_data$V1
# combat_expr_data = combat_expr_data[,-c(1)]
# colnames_expr = colnames(combat_expr_data)




t1 = read.csv('survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/R_workspace/survival/expr_data.csv',row.names = 1)

t2 = as.matrix(t2)
hist(log2(t2),breaks=100,col="blue") 

batch = t1$sequencing.protocol

batch[batch=='anchordx'] = 1
batch[batch=='novogene'] = 2
batch = t1$sequencing.protocol

batch[batch=='anchordx' & t1$gender==1] = 1
batch[batch=='anchordx' & t1$gender==2] = 2
batch[batch=='novogene' & t1$gender==1] = 3
batch[batch=='novogene' & t1$gender==2] = 4


modcombat = model.matrix(~1, data = t1)

combat_data = ComBat(dat=t2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)


aft_data = combat_data
aft_data[aft_data==0]=NA
mean_median = mean(apply(aft_data,2,median,na.rm = TRUE),na.rm = TRUE)

scale_median_3 <- function(df){
  a = median(df,na.rm = TRUE)
  c = (df-a+mean_median)
  return(c)
}

scale_median_data =  apply(combat_data_protocol,2,scale_median_3)
#a11 = apply(scale_median_data,2,median,na.rm = TRUE)




expr_data = t2

par(mfrow=c(1,2),cell=2)
hist(log2(expr_data),breaks=150,col="blue",main="Original data") 
hist(log2(combat_data),breaks=200,col="blue",main = "Remove the batch effects") 
