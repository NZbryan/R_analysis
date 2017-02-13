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

original_methy_data = fread('D:/R_workspace/survival/data_expr/广西450K甲基化matrix betas_1.csv')
patient_info = read.csv('D:/R_workspace/survival/PPT/patient_info.csv',row.names = 1)

methy_data = unite(original_methy_data, "index_met", c(2,3,4))
index_methy_data = methy_data$index_met
methy_data = methy_data[,-c(1,2,3)]

col_me = colnames(methy_data)
length_me = length(col_me)
for (i in 1:length_me)
{
  if(substr(col_me[i], 1, 1)=="2"){col_me[i] =paste0("X",col_me[i])}
}

colnames(methy_data) =col_me

# "A7","GX009","X2C2","X2D1","A8"
a11 = rownames(patient_info)
a22 = colnames(methy_data)
out_a = intersect(a11,a22)

me_patient_info = patient_info[out_a,]
af_methy_data = subset(methy_data,select = out_a)
me_patient_info$SurvObj <- with(me_patient_info, Surv(survival.month, cencor.status == 1))

pvalue_length = dim(af_methy_data)[1]
pvalue_data = array(1:pvalue_length)
hr_data = array(1:pvalue_length)

af_methy_data = t(af_methy_data)
af_methy_data = scale(af_methy_data)
af_methy_data = as.data.frame(af_methy_data)


for (i in 1:pvalue_length)
{
  model = coxph(SurvObj~af_methy_data[,i],data = me_patient_info)
  pra = glance(model)$p.value.log
  pvalue_data[i] = pra
  hr_data[i] = model$coefficients
}

hr_data = exp(hr_data)
fdr_data = p.adjust(pvalue_data, method="fdr")
pvalue_hr_fdr = data.frame(pvalue_data,hr_data,fdr_data)
rownames(pvalue_hr_fdr) = index_methy_data


# write.csv(pvalue_hr_fdr,"D:/R_workspace/survival/PPT/PPT/methy_scale.csv")




