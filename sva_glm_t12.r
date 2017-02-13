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

combat_expr_data = fread('D:/R_workspace/survival/PPT/combat_expr_data.csv')
patient_info = read.csv('D:/R_workspace/survival/PPT/patient_info.csv',row.names = 1)


combat_expr_data = as.data.frame(combat_expr_data)
rownames(combat_expr_data) = combat_expr_data$V1
combat_expr_data = combat_expr_data[,-c(1)]
colnames_expr = colnames(combat_expr_data)

patient_info$SurvObj <- with(patient_info, Surv(survival.month, cencor.status == 1))

pvalue_length = dim(combat_expr_data)[2]
pvalue_data = array(1:pvalue_length)
hr_data = array(1:pvalue_length)

combat_expr_data = scale(combat_expr_data)

for (i in 1:pvalue_length)
{
  model = coxph(SurvObj~combat_expr_data[,i],data = patient_info)
  pra = glance(model)$p.value.log
  pvalue_data[i] = pra
  hr_data[i] = model$coefficients
}



# pvalue_logrank=as.data.frame(pvalue_data)
# hr_logrank = as.data.frame(hr_data)
hr_data = exp(hr_data)
pvalue_hr = data.frame(pvalue_data,hr_data)
rownames(pvalue_hr) = colnames_expr
##write.csv(pvalue_hr,'D:/R_workspace/survival/PPT/pvalue_hr.csv')
# subset(pvalue_hr,hr_data>b11 & pvalue_data<0.05)

pvalue_logrank_005 = subset(pvalue_logrank,pvalue_data<0.005)
pvalue_logrank_005_names = rownames(pvalue_logrank_005)

x = subset(combat_expr_data,select = pvalue_logrank_005_names)
x$gender = patient_info$gender
x = as.matrix(x)
y = patient_info$survival.month



fit_cox = glmnet(x, patient_info$SurvObj, family="cox", alpha=0.1,lambda = 0.4)
pred = predict(fit_cox,x)

k11 = predict(fit_cox,x,type = "response")
k22 = data.frame(patient_info$SurvObj,k11)






set.seed(1235810)
par(mfrow=c(1,3))
cvfit_1 = cv.glmnet(x, y, family = "cox",maxit=15000,alpha=1,nlambda = 150)

cvfit_2 = cv.glmnet(x, y, family = "cox",maxit=15000,alpha=0,nlambda = 150)

cvfit_3 = cv.glmnet(x, y, family = "cox",maxit=15000,alpha=0.4,nlambda = 150)

plot(cvfit_1)
plot(cvfit_2)
plot(cvfit_3)


cvfit_3 = cv.glmnet(x, y, family = "cox",maxit=15000,alpha=0.0005,nlambda = 150)
plot(cvfit_3)


cvfit_cox = cv.glmnet(x, patient_info$SurvObj, family = "cox",maxit=100000,alpha=0.01,nlambda = 150)
plot(cvfit_cox)

fit_cox = glmnet(x, patient_info$SurvObj, family="cox", alpha=0.1,lambda = 0.4)
pred = predict(fit_cox,x)

k11 = predict(fit_cox,x,type = "response")
k22 = data.frame(patient_info$SurvObj,k11)




cvfit_cox_t1 = cv.glmnet(as.matrix(combat_expr_data), patient_info$SurvObj, family = "cox",maxit=10000,alpha=0.1,nlambda = 150)
plot(cvfit_cox_t1)



km_a = npsurv(SurvObj ~ pred, data = patient_info)
survplot(km_a,xlab = 'months', n.risk   = TRUE,cex.n.risk = 1,type="kaplan-meier",lty=1,lwd=2,
         label.curves = list(method = "arrow", cex = 1.2),col=c("orange","blue"),abbrev.label=TRUE)


cvfit_cox = cv.glmnet(x, patient_info$SurvObj, family = "cox",maxit=10000,alpha=,nlambda = 150)
plot(cvfit_cox)


for(k in 100:1000)
{
  cvfit_test = cv.glmnet(x, patient_info$SurvObj, family = "cox",maxit=10000,alpha=1/k,nlambda = 150)
  fit_test = glmnet(x, patient_info$SurvObj, family="cox", alpha=1/k,lambda = cvfit_test$lambda.min)
  if(fit_test$dev.ratio>0.7){break}
}





