require(sva)
require(ggplot2)
require(glmnet)
require(survival)
#require(reshape2)
require(rms)
require(caret)
#require(affy)
#require(affydata)

t1 = read.csv('survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/R_workspace/survival/expr_data.csv',row.names = 1)
rownames(t1) = colnames(t2)
#outlier_samples =c(A7,GX009,X2C2,X2D1)
t1 = t1[!rownames(t1) %in% c("A7","GX009","X2C2","X2D1"),]
t2 = subset(t2,select = -c(A7,GX009,X2C2,X2D1))


rm_zero_func<-function(df,...)
{
  return(all(df==0))
}
col_t2 = apply(t2,1,rm_zero_func)
t2 = t2[!col_t2,]

t2 = as.matrix(t2)
t1$gender[which(t1$gender=='F')]=1
t1$gender[which(t1$gender=='M')]=2
t1$gender = as.integer(t1$gender)

# sequencing.protocol  batch
batch = t1$sequencing.protocol

batch[batch=='anchordx' & t1$gender==1] = 1
batch[batch=='anchordx' & t1$gender==2] = 2
batch[batch=='novogene' & t1$gender==1] = 3
batch[batch=='novogene' & t1$gender==2] = 4

modcombat = model.matrix(~1, data = t1)
combat_data_protocol  = ComBat(dat=t2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

aft_data = combat_data_protocol
aft_data[aft_data==0]=NA
mean_median = mean(apply(aft_data,2,median,na.rm = TRUE),na.rm = TRUE)

scale_median_3 <- function(df){
  a = median(df,na.rm = TRUE)
  c = (df-a+mean_median)
  return(c)
}

scale_median_data =  apply(combat_data_protocol,2,scale_median_3)
cv_data = t(scale_median_data)

## pvalues < 0.005
pValuesComBat_func = function(df_adjust,df_expr,...)
{
  mod = model.matrix(~as.factor(cencor.status), data=df_adjust)
  mod0 = model.matrix(~1, data = df_adjust)
  
  df_expr = t(df_expr)
  pValuesComBat = f.pvalue(df_expr,mod,mod0)
  
  qValuesComBat = p.adjust(pValuesComBat,method="BH")
  
  pValuesComBat_005 = pValuesComBat[pValuesComBat<0.05]
  pValuesComBat_005 = rownames(data.frame(pValuesComBat_005))
  return(pValuesComBat_005)
}




## devide data to three parts
set.seed(201711600)
cv_Folds = createFolds(1:dim(cv_data)[1],k=5)
#png("D:/R_workspace/survival/PPT/cv5.png",width = 2560,height = 1440)
par(mfrow = c(5,1))

## for
for(k in 1:length(cv_Folds))
{
  data_adjust_1 = t1[cv_Folds[[k]],]
  data_adjust_others = t1[-cv_Folds[[k]],]
  
  data_expr_1 = cv_data[cv_Folds[[k]],]
  data_expr_others = cv_data[-cv_Folds[[k]],]
  
  sample_data = as.data.frame(data_expr_others)
  
  pValuesComBat_005 = pValuesComBat_func(data_adjust_others,data_expr_others)
  pv = sample_data[pValuesComBat_005]
  pv$gender = data_adjust_others$gender
  
  ## set the matircs of x and y
  x = as.matrix(pv)
  y = data_adjust_others$survival.month
  
  ## train the model
  fit_model = glmnet(x, y, family="gaussian", alpha=0.42,lambda = 1.3750)
  
  ## choose the features
  a11 = fit_model$beta
  beta = a11[which(a11 != 0)]
  chosen_names = rownames(a11)[which(a11 != 0)]
  chosen_features = data.frame(chosen_names,beta)
  
  data_chosen_m = as.data.frame(x)
  data_chosen = data_chosen_m[chosen_names]
  data_chosen_x = as.matrix(data_chosen)
  
  ## predict
  data_expr_1 = as.data.frame(data_expr_1)
  data_expr_1 = data_expr_1[pValuesComBat_005]
  data_expr_1$gender = data_adjust_1$gender

  
  testData_x = as.matrix(data_expr_1)
  testData_y = data_adjust_1$survival.month
  
  y_predict = predict(fit_model,testData_x)

  data_a = data.frame(testData_y,y_predict)
  colnames(data_a)=c('y_true','y_predict')
  data_a$y_true_status = data_adjust_1$cencor.status
  
  # log rank & km curve
  data_a$predict_binary = data_a$y_predict
  md_va = data_a$predict_binary>39
  data_a$predict_binary[md_va] = 0
  data_a$predict_binary[!md_va] = 1
  status_mid1 = table(data_a$predict_binary,data_a$y_true_status)[1]+
    table(data_a$predict_binary,data_a$y_true_status)[4]
  status_mid2 = table(data_a$predict_binary,data_a$y_true_status)[2]+
    table(data_a$predict_binary,data_a$y_true_status)[3]
  if(status_mid1<status_mid2)
  {
    data_a$predict_binary[data_a$predict_binary == 1] = 2
    data_a$predict_binary[data_a$predict_binary == 0] = 1
    data_a$predict_binary[data_a$predict_binary == 2] = 0
  }

  SurvObj_a <- with(data_a, Surv(y_true, y_true_status == 1))
  km_a = npsurv(SurvObj_a ~ predict_binary, data = data_a)
  
  log_rank_predict_a=survdiff(SurvObj_a ~ predict_binary, data = data_a)
  log_rank_predict_a_pvalue = pchisq(log_rank_predict_a$chisq, length(log_rank_predict_a$n)-1, lower.tail = FALSE)
  
  survplot(km_a,xlab = 'months', n.risk   = TRUE,cex.n.risk = 1,type="kaplan-meier",lty=1,lwd=2,
           label.curves = list(method = "arrow", cex = 1.2),col=c("orange","blue"),abbrev.label=TRUE)
  text(60,1,sprintf("p-value of log-rank test: %.10f", log_rank_predict_a_pvalue),col='brown2')
  title_K = paste("test",k)
  title(title_K)
  
}

#dev.off()




