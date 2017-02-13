require(sva)
require(data.table)
require(ggplot2)
require(glmnet)
require(survival)
#require(reshape2)
require(rms)
require(caret)
require(broom)
require(stringr)
require(tidyr)
#
original_methy_data = fread('D:/R_workspace/survival/PPT/methy_data/60例样本甲基化芯片beta矩阵betas_1.csv')
patient_info = read.csv('D:/R_workspace/survival/PPT/patient_info.csv',row.names = 1,stringsAsFactors = FALSE)

methy_data = unite(original_methy_data, "index_met", c(2,3,4))
index_methy_data = methy_data$index_met
methy_data = methy_data[,-c(1,2,3)]

col_me = colnames(methy_data)
length_me = length(col_me)
for (i in 1:length_me)
{
  if(substr(col_me[i], 1, 1)=="2"){col_me[i] =paste0("X",col_me[i])}
  if(substr(col_me[i], 1, 5)=="GX001"){col_me[i] = paste0("GX",str_split_fixed(col_me[i],"-",2)[2])}
}
colnames(methy_data) =col_me

# "A7","GX009","X2C2","X2D1","A8"
a11 = rownames(patient_info)
a22 = colnames(methy_data)
out_a = intersect(a11,a22)

me_patient_info = patient_info[out_a,]
af_methy_data = subset(methy_data,select = out_a)
af_methy_data = as.data.frame(af_methy_data)
rownames(af_methy_data) = index_methy_data
af_methy_data = na.omit(af_methy_data)
# write.csv(me_patient_info,"D:/R_workspace/survival/PPT/methy_data/60data/60me_patient_info.csv")

me_patient_info$SurvObj <- with(me_patient_info, Surv(survival.month, cencor.status == 1))

# sequencing.protocol  batch
batch = me_patient_info$sequencing.protocol
batch[batch=='anchordx'] = 1
batch[batch=='novogene'] = 2
batch = as.integer(batch)

modcombat = model.matrix(~1, data = me_patient_info)
af_methy_data = as.matrix(af_methy_data)

combat_data_protocol  = ComBat(dat=af_methy_data, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)
# write.csv(combat_data_protocol,"D:/R_workspace/survival/PPT/methy_data/60data/60methy_combat_data_protocol.csv")

# coxph for single gene
pvalue_length = dim(combat_data_protocol)[1]
pvalue_data = array(1:pvalue_length)
hr_data = array(1:pvalue_length)

combat_data_protocol = t(combat_data_protocol)
combat_data_protocol = scale(combat_data_protocol)
combat_data_protocol = as.data.frame(combat_data_protocol)
index_combat_data = colnames(combat_data_protocol)
#
for (i in 1:pvalue_length)
{
  model = coxph(SurvObj~combat_data_protocol[,i],data = me_patient_info)
  pra = glance(model)$p.value.log
  pvalue_data[i] = pra
  hr_data[i] = model$coefficients
}

hr_data = exp(hr_data)
fdr_data = p.adjust(pvalue_data, method="fdr")
pvalue_hr_fdr = data.frame(pvalue_data,hr_data,fdr_data)
rownames(pvalue_hr_fdr) = index_combat_data

#write.csv(pvalue_hr_fdr,"D:/R_workspace/survival/PPT/methy_data/60data/combat_methy_pvalue_hr_fdr.csv")