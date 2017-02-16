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
set.seed(2017215)
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

patient_info$SurvObj <- with(patient_info, Surv(survival.month, cencor.status == 1))
# methy_combat_expr_data = na.omit(methy_combat_expr_data)

x3 = t(methy_combat_expr_data)
x3 = scale(x3)
y = patient_info$SurvObj
#Parallel computing & cv.glmnet
require(doParallel)
registerDoParallel(cores=8)
cvfit_cox3 = cv.glmnet(x3, y, family = "cox",alpha=0.15,nlambda = 150,standardize = F,parallel=TRUE)
stopImplicitCluster()
## train the model
fit_cox3 = glmnet(x3,y, family="cox",alpha=0.15,lambda = cvfit_cox3$lambda.min ,standardize = F)

a11 = fit_cox3$beta
a22 = as.matrix(a11)
a33 = as.data.frame(a22)
chosen_features = subset(a33,s0!=0)
chosen_features$index_methy = rownames(chosen_features)
chosen_features_phf_all = pvalue_hr[chosen_features$index_methy,]
chosen_features_phf_all$beta = chosen_features$s0
# write.csv(chosen_features_phf_all,"R/R_data/output/chosen_features_phf_all.csv")
# 0.05 -----------------------------------------------------------------------------


pvalue_hr_05 = subset(pvalue_hr,pvalue_data<0.05)
pvalue_hr_05_names = rownames(pvalue_hr_05)
methy_combat_expr_data = methy_combat_expr_data[pvalue_hr_05_names,]

x1 = t(methy_combat_expr_data)
x1 = scale(x1)
#Parallel computing & cv.glmnet
require(doParallel)
registerDoParallel(cores=8)
cvfit_cox1 = cv.glmnet(x1, y, family = "cox",alpha=0.15,nlambda = 150,standardize = F,parallel=TRUE)
stopImplicitCluster()
## train the model
fit_cox1 = glmnet(x1,y, family="cox",alpha=0.15,lambda = cvfit_cox1$lambda.min ,standardize = F)

a11 = fit_cox1$beta
a22 = as.matrix(a11)
a33 = as.data.frame(a22)
chosen_features = subset(a33,s0!=0)
chosen_features$index_methy = rownames(chosen_features)
chosen_features_phf_005 = pvalue_hr[chosen_features$index_methy,]
chosen_features_phf_005$beta = chosen_features$s0

# 0.005 -----------------------------------------------------------------------------
methy_combat_expr_data = fread('R/R_data/60methy_combat_data_protocol.csv')
methy_combat_expr_data = as.data.frame(methy_combat_expr_data)
rownames(methy_combat_expr_data) = methy_combat_expr_data$V1
methy_combat_expr_data = methy_combat_expr_data[,-c(1)]

pvalue_hr_005 = subset(pvalue_hr,pvalue_data<0.005)
pvalue_hr_005_names = rownames(pvalue_hr_005)
methy_combat_expr_data = methy_combat_expr_data[pvalue_hr_005_names,]

x2 = t(methy_combat_expr_data)
x2 = scale(x2)
#Parallel computing & cv.glmnet
require(doParallel)
registerDoParallel(cores=8)
cvfit_cox2 = cv.glmnet(x2, y, family = "cox",alpha=0.15,nlambda = 150,standardize = F,parallel=TRUE)
stopImplicitCluster()
## train the model
fit_cox2 = glmnet(x2,y, family="cox",alpha=0.15,lambda = cvfit_cox2$lambda.min ,standardize = F)

a11 = fit_cox2$beta
a22 = as.matrix(a11)
a33 = as.data.frame(a22)
chosen_features = subset(a33,s0!=0)
chosen_features$index_methy = rownames(chosen_features)
chosen_features_phf_0000005 = pvalue_hr[chosen_features$index_methy,]
chosen_features_phf_0000005$beta = chosen_features$s0
# write.csv(chosen_features_phf_0000005,"R/R_data/output/chosen_features_phf_005.csv")

#  common methylreads ------------------------------------------------------------------------------
chosen_features_phf_005_name = rownames(chosen_features_phf_005)
chosen_features_phf_all_name = rownames(chosen_features_phf_all)

co_name = intersect(chosen_features_phf_005_name,chosen_features_phf_all_name)


chosen_features_phf_0000005_name = rownames(chosen_features_phf_0000005)


co_name_12 = intersect(chosen_features_phf_all_name,chosen_features_phf_005_name)
co_name_13 = intersect(chosen_features_phf_all_name,chosen_features_phf_0000005_name)
co_name_23 = intersect(chosen_features_phf_005_name,chosen_features_phf_0000005_name)

asas = c("chr17_73775874_73775875","chr2_61372256_61372257","chr6_126661238_126661239","chr6_5261091_5261092")

grep("chr6_5261091_5261092",chosen_features_phf_0000005_name)


