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
#require(affy)
#require(affydata)
# options(warn=-1) 
methy_combat_expr_data = fread('D:/R_workspace/survival/PPT/methy_data/60data/60methy_combat_data_protocol.csv')
patient_info = read.csv('D:/R_workspace/survival/PPT/methy_data/60data/60me_patient_info.csv',row.names = 1,stringsAsFactors = FALSE)
# pvalue_hr = fread('D:/R_workspace/survival/PPT/methy_data/combat_methy_pvalue_hr_fdr.csv')
pvalue_hr = read.csv('D:/R_workspace/survival/PPT/methy_data/60data/combat_methy_pvalue_hr_fdr.csv',row.names = 1,stringsAsFactors = FALSE)

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
cvfit_cox1 = cv.glmnet(x, y, family = "cox",alpha=1,nlambda = 150)
cvfit_cox2 = cv.glmnet(x, y, family = "cox",alpha=0,nlambda = 150)
cvfit_cox3 = cv.glmnet(x, y, family = "cox",alpha=0.5,nlambda = 150)

par(mfrow=c(3,1))
plot(cvfit_cox1)
plot(cvfit_cox2)
plot(cvfit_cox3)



cl <- makeCluster(6)
registerDoParallel(cl)
cvfit_cox4 = cv.glmnet(x, y, family = "cox",alpha=0.4,nlambda = 150,parallel=TRUE)
stopCluster(cl)

cvfit_cox5 = cv.glmnet(x, y, family = "cox",alpha=0.6,nlambda = 150)


cvfit_cox3$glmnet.fit


plot(cvfit_cox4, xlim = c(-1.1, -0.5), ylim = c(7, 11))
