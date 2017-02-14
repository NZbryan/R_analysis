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


#Parallel computing & cv.glmnet
require(doParallel)
registerDoParallel(cores=8)
cvfit_cox4 = cv.glmnet(x, y, family = "cox",alpha=0.15,nlambda = 150,standardize = F,parallel=TRUE)
stopImplicitCluster()
## train the model
fit_cox = glmnet(x,y, family="cox",alpha=0.15,lambda = cvfit_cox4$lambda.min,standardize = F)
# coef(fit_cox_5,s=c(fit_cox_5$lambda))
# tidy(cvfit_cox4)[tidy(cvfit_cox4)$lambda==cvfit_cox4$lambda.min,]$estimate
a11 = fit_cox$beta
a22 = as.matrix(a11)
a33 = as.data.frame(a22)
chosen_features = subset(a33,s0!=0)


y_predict = predict(fit_cox,x,type = "response")

data_a = data.frame(y_predict,patient_info$SurvObj)

data_a$signature = data_a$s0
md_va = data_a$signature>1.5
data_a$signature[md_va] = "High-risk"
data_a$signature[!md_va] = "Low-risk"

km_a = npsurv(patient_info.SurvObj ~ signature, data = data_a)

log_rank_predict_a=survdiff(patient_info.SurvObj ~ signature, data = data_a)
log_rank_predict_a_pvalue = pchisq(log_rank_predict_a$chisq, length(log_rank_predict_a$n)-1, lower.tail = FALSE)

par(mar=c(5.1,4.1,4.1,7.1),cex=1.2)
survplot(km_a,xlab = 'months', n.risk   = TRUE,y.n.risk=0,cex.n.risk = 1,type="kaplan-meier",lty=1,lwd=3,
         label.curves = list(method = "arrow", cex = 1.2),col=c("orange","blue"),abbrev.label=FALSE)
text(60,0.65,sprintf("p-value of log-rank test: %.10f", log_rank_predict_a_pvalue),col='brown2')
title_K = paste("test:",1)
title(title_K)
# a11 = fit_cox$beta
# beta = as.matrix(a11)
# beta = as.data.frame(beta)
# chosen_names = rownames(beta)
# 
# chosed_p_hr = pvalue_hr[pvalue_logrank_005_names,]
# add_gender = c(0.0291,2.9888)
# chosed_p_hr = rbind(chosed_p_hr,add_gender)
# rownames(chosed_p_hr)[421] = "gender"
# chosed_p_hr_beta = chosed_p_hr
# chosed_p_hr_beta$beta = beta$s0
# chosed_p_hr_beta = chosed_p_hr_beta[order(-beta),]
# 



