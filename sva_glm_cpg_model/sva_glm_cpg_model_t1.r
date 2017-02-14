require(sva)
require(data.table)
require(ggplot2)
require(glmnet)
require(survival)
#require(reshape2)
require(rms)
require(caret)
require(broom)
require(tidyr)
require(openxlsx)
require(stringr)
# options(warn=-1) 
# rm(list=ls())
methy_combat_expr_data = fread('R/R_data/60methy_combat_data_protocol.csv')
patient_info = read.csv('R/R_data/60me_patient_info.csv',row.names = 1,stringsAsFactors = FALSE)
# pvalue_hr = fread('D:/R_workspace/survival/PPT/methy_data/combat_methy_pvalue_hr_fdr.csv')
pvalue_hr = read.csv('R/R_data/combat_methy_pvalue_hr_fdr.csv',row.names = 1,stringsAsFactors = FALSE)

methy_combat_expr_data = as.data.frame(methy_combat_expr_data)
rownames(methy_combat_expr_data) = methy_combat_expr_data$V1
methy_combat_expr_data = methy_combat_expr_data[,-c(1)]
# colnames_expr = colnames(methy_combat_expr_data)

pvalue_hr_05 = subset(pvalue_hr,pvalue_data<0.05)
pvalue_hr_05_names = rownames(pvalue_hr_05)

methy_combat_expr_data = methy_combat_expr_data[pvalue_hr_05_names,]
# x,y
patient_info$SurvObj <- with(patient_info, Surv(survival.month, cencor.status == 1))
# methy_combat_expr_data = na.omit(methy_combat_expr_data)

x = t(methy_combat_expr_data)
x = scale(x)
y = patient_info$SurvObj




# model
cvfit_cox1 = cv.glmnet(x, y, family = "cox",alpha=1,nlambda = 150,,standardize = F)
cvfit_cox2 = cv.glmnet(x, y, family = "cox",alpha=0,nlambda = 150,standardize = F)
cvfit_cox3 = cv.glmnet(x, y, family = "cox",alpha=0.5,nlambda = 150,standardize = F)

par(mfrow=c(3,1))
plot(cvfit_cox1)
plot(cvfit_cox2)
plot(cvfit_cox3)



cvfit_cox4 = cv.glmnet(x, y, family = "cox",alpha=0.15,nlambda = 150)



cvfit_cox = cv.glmnet(x, y, family = "cox",alpha=0.15,nlambda = 150)
#plot(cvfit_cox)

## train the model
fit_cox = glmnet(x,y, family="cox",alpha=0.15,lambda = cvfit_cox4$lambda.min,standardize = F)



a11 = fit_cox$beta
beta = as.matrix(a11)
beta = as.data.frame(beta)
chosen_names = rownames(beta)

chosed_p_hr = pvalue_hr[pvalue_logrank_005_names,]
add_gender = c(0.0291,2.9888)
chosed_p_hr = rbind(chosed_p_hr,add_gender)
rownames(chosed_p_hr)[421] = "gender"
chosed_p_hr_beta = chosed_p_hr
chosed_p_hr_beta$beta = beta$s0
chosed_p_hr_beta = chosed_p_hr_beta[order(-beta),]



fit_cox_5 = glmnet(x,y, family="cox",alpha=0.15,nlambda = 150,standardize = F)

coef(fit_cox_5,s=c(fit_cox_5$lambda))

#Parallel computing
require(doParallel)
registerDoParallel(cores=8)
cvfit_cox4 = cv.glmnet(x, y, family = "cox",alpha=0.15,nlambda = 150,standardize = F,parallel=TRUE)
stopImplicitCluster()


tidy(cvfit_cox4)[tidy(cvfit_cox4)$lambda==cvfit_cox4$lambda.min,]$estimate


## alpha = i/10

require(doParallel)
registerDoParallel(cores=8)
for(i in 1:9)
{
  cvfit_test = cv.glmnet(x, y, family = "cox",alpha=i/10,nlambda = 150,standardize = F,parallel=TRUE)
  output = tidy(cvfit_test)[tidy(cvfit_test)$lambda==cvfit_test$lambda.min,]$estimate
  print(output)
}
stopImplicitCluster()

###output ：  estimate (median) of mean-squared error or other criterion

# [1] 8.475498
# [1] 8.806469
# [1] 8.682234
# [1] 8.508602
# [1] 8.988216
# [1] 9.195006
# [1] 8.486525
# [1] 8.895532
# [1] 8.638769


## alpha = i/20

require(doParallel)
registerDoParallel(cores=8)
for(i in 1:9)
{
  cvfit_test = cv.glmnet(x, y, family = "cox",alpha=i/20,nlambda = 150,standardize = F,parallel=TRUE)
  output = tidy(cvfit_test)[tidy(cvfit_test)$lambda==cvfit_test$lambda.min,]$estimate
  print(output)
}
stopImplicitCluster()

###output ：

# [1] 8.659479
# [1] 8.563464
# [1] 8.246525
# [1] 8.874835
# [1] 8.696664
# [1] 9.010895
# [1] 8.717095
# [1] 8.778384
# [1] 8.919419
