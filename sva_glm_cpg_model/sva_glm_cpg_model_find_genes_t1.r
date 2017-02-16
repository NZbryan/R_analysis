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
set.seed(2017215)
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
chosen_features$index_methy = rownames(chosen_features)
chosen_features_phf = pvalue_hr[chosen_features$index_methy,]
chosen_features_phf$beta = chosen_features$s0
# write.csv(chosen_features_phf,"R/R_data/output/chosen_features_phf005.csv")

### input cpg_to_genes table
genes_data2 = read.csv('R/R_data/multianno.txt',stringsAsFactors = FALSE,sep="")
gene_t1 = genes_data2[c(1,2,3,6,7)]
gene_t2 = unite(gene_t1, "index_methy", c(1,2,3))
gene_t2$count_nchar = apply(gene_t2["Gene.ensGene"],1,nchar)
gene_t3 = subset(gene_t2,count_nchar == 15)

co_cpg = intersect(chosen_features$index_methy,gene_t3$index_methy)
gene_t4 = gene_t3[gene_t3$index_methy %in% co_cpg,]

co_data = merge(x = gene_t4, y = chosen_features, by = "index_methy", all.x=TRUE)
co_data = subset(co_data,select = -count_nchar)

chosed_pvalue_hr = pvalue_hr
chosed_pvalue_hr$index_methy = rownames(chosed_pvalue_hr)
colnames(chosed_pvalue_hr)[1:3] = c("p-value of cpg","hr of cpg","fdr of cpg")

co_data_addphf = merge(x = co_data, y = chosed_pvalue_hr, by = "index_methy", all.x=TRUE)
### input genes model features

gene_model_features = read.xlsx('R/R_data/chosed_p_hr_beta_new.xlsx')
colnames(gene_model_features) = c("features","p-value of RNASep","hr of RNASep","fdr of RNASep")
gene_model_features = gene_model_features[-c(3),]

hclust_name = function(df)
{ 
  s1 = str_split_fixed(df,"_",2)[1]
}

Gene.ensGene = apply(gene_model_features["features"],1,hclust_name)

gene_model_features_t2 = cbind(gene_model_features[,1],Gene.ensGene,gene_model_features[,2:ncol(gene_model_features)])

## co
co_Gene.ensGene = intersect(co_data_addphf$Gene.ensGene,as.character(gene_model_features_t2$Gene.ensGene))


gene_model_features_t3 = gene_model_features_t2[gene_model_features_t2$Gene.ensGene %in% co_Gene.ensGene,]
# gene_model_features[, 1]    Gene.ensGene p-value of RNASep hr of RNASep fdr of RNASep
# 23     ENSG00000132475_H3F3B ENSG00000132475      0.0048562714    1.6238568    0.05246448
# 49     ENSG00000214113_LYRM4 ENSG00000214113      0.0017911062    1.8212266    0.01507438
# 50     ENSG00000203760_CENPW ENSG00000203760      0.0005033875    1.8010359    0.01467221
# 104  ENSG00000237651_C2orf74 ENSG00000237651      0.0017620926    0.5254502   -0.05865387


co_data_addphf_t2 = co_data_addphf[co_data_addphf$Gene.ensGene %in% co_Gene.ensGene,]
# index_methy Func.ensGene    Gene.ensGene           s0 p-value of cpg hr of cpg fdr of cpg
# 91   chr17_73775874_73775875     intronic ENSG00000132475 -0.006937409   3.264957e-04 0.4041507 0.02593493
# 139   chr2_61372256_61372257         UTR5 ENSG00000237651  0.003928449   4.694580e-04 2.2288583 0.02879789
# 169 chr6_126661238_126661239     upstream ENSG00000203760 -0.001440588   1.419661e-03 0.4225360 0.04117087
# 182     chr6_5261091_5261092         UTR5 ENSG00000214113 -0.030895643   1.991186e-05 0.2983138 0.01509267

### just 4 common genes between RNASeq and cpg
  