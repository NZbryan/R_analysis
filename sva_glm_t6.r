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
#a11 = apply(scale_median_data,2,median,na.rm = TRUE)


mod = model.matrix(~as.factor(cencor.status), data=t1)
mod0 = model.matrix(~1, data = t1)

pValuesComBat = f.pvalue(scale_median_data,mod,mod0)

qValuesComBat = p.adjust(pValuesComBat,method="BH")

pValuesComBat_005 = pValuesComBat[pValuesComBat<0.05]
pValuesComBat_005 = rownames(data.frame(pValuesComBat_005))

sample_data = as.data.frame(t(scale_median_data))

pv = sample_data[pValuesComBat_005]
pv$gender = t1$gender

x = as.matrix(pv)
y = t1$survival.month





#set.seed(12356896)
#cvfit_4 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0.42,nlambda = 150)
#plot(cvfit_4)
#chosen_index = rownames(coefficients)[which(coefficients != 0)]

fit_04 = glmnet(x, y, family="gaussian", alpha=0.42,lambda = 1.3750)

a11 = fit_04$beta


beta = a11[which(a11 != 0)]
chosen_names = rownames(a11)[which(a11 != 0)]
chosen_features = data.frame(chosen_names,beta)
#write.csv(chosen_features,'D:/R_workspace/survival/PPT/chosen_features_1.csv')



data_chosen_m = as.data.frame(x)
data_chosen = data_chosen_m[chosen_names]
data_chosen_x = as.matrix(data_chosen)
#fit_05 = glmnet(data_chosen_x, y, family="gaussian", alpha=0.42,lambda = 1.3750)

pdf('D:/R_workspace/survival/PPT/test1.pdf',width=15,height = 12)
for(k in 1:10){
  num_set = as.integer(paste(171120,k,sep=''))
  set.seed(num_set)
  cv_Folds = createFolds(1:76,k=2)
  
  trainData_x = data_chosen_x[cv_Folds$Fold1,]
  testData_x = data_chosen_x[cv_Folds$Fold2,]
  trainData_y = y[cv_Folds$Fold1]
  testData_y = y[cv_Folds$Fold2]
  
  
  ## a组
  fit_a = glmnet(trainData_x, trainData_y, family="gaussian", alpha=0.42,lambda = 1.3750)
  pre_y_a = predict(fit_a,testData_x)
  
  data_a = data.frame(testData_y,pre_y_a)
  colnames(data_a)=c('y_a','predict_a')
  
  y_a_status = t1$cencor.status[cv_Folds$Fold2]
  data_a$y_a_status = y_a_status
  
  ## b组
  fit_b = glmnet(testData_x, testData_y, family="gaussian", alpha=0.42,lambda = 1.3750)
  pre_y_b = predict(fit_b,trainData_x)
  
  data_b = data.frame(trainData_y,pre_y_b)
  colnames(data_b)=c('y_b','predict_b')
  
  y_b_status = t1$cencor.status[cv_Folds$Fold1]
  data_b$y_b_status = y_b_status
  
  # log rank & km curve
  data_a$predict_data = data_a$predict_a
  data_a$predict_a[data_a$predict_a<=39] = 1
  data_a$predict_a[data_a$predict_a>39] = 2
  
  data_b$predict_data = data_b$predict_b
  data_b$predict_b[data_b$predict_b<=39] = 1
  data_b$predict_b[data_b$predict_b>39] = 2
  
  SurvObj_a <- with(data_a, Surv(y_a, y_a_status == 1))
  km_a = npsurv(SurvObj_a ~ predict_a, data = data_a)
  
  SurvObj_b <- with(data_b, Surv(y_b, y_b_status == 1))
  km_b = npsurv(SurvObj_b ~ predict_b, data = data_b)
  
  log_rank_predict_a=survdiff(SurvObj_a ~ predict_a, data = data_a)
  log_rank_predict_a_pvalue = pchisq(log_rank_predict_a$chisq, length(log_rank_predict_a$n)-1, lower.tail = FALSE)
  
  log_rank_predict_b=survdiff(SurvObj_b ~ predict_b, data = data_b)
  log_rank_predict_b_pvalue = pchisq(log_rank_predict_b$chisq, length(log_rank_predict_b$n)-1, lower.tail = FALSE)
  title_a = paste("Group a : ",k)
  title_b = paste("Group b : ",k)
  
  par(mfrow=c(2,1))
  survplot(km_a,xlab = 'months', n.risk   = TRUE,cex.n.risk = 1,type="kaplan-meier",
           label.curves = list(method = "arrow", cex = 1.2),col=c(1:4),abbrev.label=TRUE)
  text(90,0.8,sprintf("p-value of log-rank test: %.10f", log_rank_predict_a_pvalue),col='brown2')
  title(title_a)
  
  survplot(km_b,xlab = 'months', n.risk   = TRUE,cex.n.risk = 1,type="kaplan-meier",
           label.curves = list(method = "arrow", cex = 1.2),col=c(1:4),abbrev.label=TRUE)
  text(90,0.8,sprintf("p-value of log-rank test: %.10f", log_rank_predict_b_pvalue),col='brown2')
  title(title_b)
}

dev.off()

## 画图





num_set = 1711202

for(k in 1:10){
  num_set = as.integer(paste(171120,k,sep=''))
  set.seed(num_set)
  cv_Folds = createFolds(1:76,k=2)
  print(cv_Folds)
}




#par(mfrow=c(2,1))
plot_grid(plotlist = plots)
p = ggplot(data = data_a, aes(x = y_a, y = predict_data))+ geom_point(col='blue')
p + geom_abline(intercept = 0, slope = 1,col = 'red',size=1.5)+geom_smooth(method=lm)

p = ggplot(data = data_b, aes(x = y_b, y = predict_data))+ geom_point(col='blue')
p + geom_abline(intercept = 0, slope = 1,col = 'red',size=1.5)+
  geom_smooth(method=lm)








