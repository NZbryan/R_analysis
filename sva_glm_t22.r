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
##write.csv(pvalue_hr,'D:/R_workspace/survival/PPT/pvalue_hr_new.csv')
# subset(pvalue_hr,hr_data>b11 & pvalue_data<0.05)

pvalue_logrank_005 = subset(pvalue_hr,pvalue_data<0.005)
pvalue_logrank_005_names = rownames(pvalue_logrank_005)


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


Significant_Predictive_Genes = chosed_p_hr_beta


