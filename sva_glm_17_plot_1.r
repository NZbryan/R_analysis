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
#

pgd = read.csv('D:/R_workspace/survival/PPT/methy_data/60data/outp/1.csv',row.names = 1,stringsAsFactors = FALSE)

xx1 = pgd["pvalue_data"]
xx1$type = "p-value"

xx2 = pgd["fdr_data"]
xx2$type = "fdr"

names(xx2)[1] = names(xx1)[1]
co_1 = rbind(xx1,xx2)
colnames(co_1)[1] = "value"

ggplot(co_1)+aes(x=value)+geom_histogram(position = 'dodge',binwidth=0.01)+aes(fill=type)



combat_expr_data = fread('D:/R_workspace/survival/PPT/combat_expr_data.csv')
patient_info = read.csv('D:/R_workspace/survival/PPT/methy_data/60data/outp/patient_info.csv',row.names = 1)
pvalue_logrank = read.csv('D:/R_workspace/survival/PPT/pvalue_logrank.csv',row.names = 1)
# pvalue_hr = read.csv('D:/R_workspace/survival/PPT/pvalue_hr.csv',row.names = 1)
pvalue_hr = read.csv('D:/R_workspace/survival/PPT/pvalue_hr_new.csv',row.names = 1)


combat_expr_data = as.data.frame(combat_expr_data)
rownames(combat_expr_data) = combat_expr_data$V1
combat_expr_data = combat_expr_data[,-c(1)]
colnames_expr = colnames(combat_expr_data)

#IFNG 、IFNG-AS1 、WNT5A
# mid1 = grep(paste0("_",y2[i],"$"),colnames_expr)

require(openxlsx)
t11 = read.xlsx("D:/R_workspace/survival/PPT/methy_data/60data/outp/t11.xlsx")
rownames(t11) = t11$SampleID

t12 = t11[c(1)]

t12_names = rownames(t12)
length_me = length(t12_names)
for (i in 1:length_me)
{
  if(substr(t12_names[i], 1, 1)=="2"){t12_names[i] =paste0("X",t12_names[i])}
}
rownames(t12) =t12_names


co_names = intersect(rownames(patient_info),t12_names)
co_data = t12[co_names,]
patient_info$Group = co_data


find_genes = array(NA,3)
gene_data = c("IFNG","IFNG-AS1","WNT5A")
for (i in 1:3)
{
  mid1 = grep(paste0("_",gene_data[i],"$"),colnames_expr)
  find_genes[i] = colnames_expr[mid1]
}

three_gene_data = combat_expr_data[find_genes]


three_gene_data$Group = co_data




plotmatrix(three_gene_data)+geom_smooth()



ENSG00000111537_IFNG = three_gene_data[c("ENSG00000111537_IFNG","Group")]
ENSG00000255733_IFNG_AS1 = three_gene_data[c("ENSG00000255733_IFNG-AS1","Group")]
ENSG00000114251_WNT5A = three_gene_data[c("ENSG00000114251_WNT5A","Group")]


a11 = ENSG00000111537_IFNG
a22 = ENSG00000255733_IFNG_AS1
a33 = ENSG00000114251_WNT5A

colnames(ENSG00000111537_IFNG)[1] = "value"
colnames(ENSG00000255733_IFNG_AS1)[1] = "value"
colnames(ENSG00000114251_WNT5A)[1] = "value"


a11$color = "ENSG00000111537_IFNG"
a22$color = "ENSG00000255733_IFNG_AS1"
a33$color = "ENSG00000114251_WNT5A"

colnames(a11)[1] = "value"
colnames(a22)[1] = "value"
colnames(a33)[1] = "value"

co_point_data = rbind(a11,a22,a33)

qplot(Group,value,data=co_point_data,geom = "point",color=color,size=0.5)


ENSG00000111537_IFNG$color = "orange"
ENSG00000255733_IFNG_AS1$color = "blue"
ENSG00000114251_WNT5A$color = "green"

p1<-qplot(Group,value,data=ENSG00000111537_IFNG,geom = "point",size=0.5,color=I("steelblue"),main="ENSG00000111537_IFNG")+
  theme(plot.title = element_text(hjust = 0.5))
p2<-qplot(Group,value,data=ENSG00000255733_IFNG_AS1,geom = "point",colour=I("orange"),size=0.5,main="ENSG00000255733_IFNG_AS1")+
  theme(plot.title = element_text(hjust = 0.5))
p3<-qplot(Group,value,data=ENSG00000114251_WNT5A,geom = "point",size=0.5,colour=I("limegreen"),main="ENSG00000114251_WNT5A")+
  theme(plot.title = element_text(hjust = 0.5))
multiplot(p1, p2, p3, cols=3)



qplot(Group,value,data=ENSG00000111537_IFNG,color=1)

fdr_data = p.adjust(pvalue_hr$pvalue_data, method="fdr")
pvalue_hr$fdr_data =fdr_data

# write.csv(pvalue_hr,"D:/R_workspace/survival/PPT/methy_data/60data/outp/pvalue_hr.csv")
