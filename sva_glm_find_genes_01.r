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

genes_data2 = read.csv('D:/R_workspace/survival/PPT/methy_data/multianno.txt',stringsAsFactors = FALSE,sep="")


HALLMARK_INTERFERON_GAMMA_RESPONSE = read.xlsx('D:/R_workspace/survival/PPT/methy_data/60data/find_gene_data/HALLMARK_INTERFERON_GAMMA_RESPONSE.xlsx')
nature09247 = read.xlsx('D:/R_workspace/survival/PPT/methy_data/60data/find_gene_data/nature09247_s4.xlsx')


feature = read.xlsx('D:/R_workspace/survival/PPT/methy_data/RNASep.xlsx')
hclust_name = function(df)
{ 
  s1 = str_split_fixed(df,"_",2)[2]
}

b2 = apply(feature,1,hclust_name)

feature$names = b2


pvalue_hr = read.csv('D:/R_workspace/survival/PPT/methy_data/60data/outp/pvalue_hr.csv',row.names = 1)


names_all_gens = rownames(pvalue_hr)
asa = str_split_fixed(names_all_gens,"_",2)

all_gens_colu = data.frame(asa)
all_gens_symbol = all_gens_colu$X2

pvalue_hr$symbol = all_gens_symbol

HALLMARK_INTERFERON_GAMMA_RESPONSE_GENE.SYMBOL = HALLMARK_INTERFERON_GAMMA_RESPONSE$GENE.SYMBOL

nature09247_GENE.SYMBOL = nature09247$Symbol

co_all_HALLMARK = intersect(HALLMARK_INTERFERON_GAMMA_RESPONSE_GENE.SYMBOL,all_gens_symbol)

co_nature09247 = intersect(nature09247_GENE.SYMBOL,all_gens_symbol)


co_all_HALLMARK_phf = pvalue_hr[pvalue_hr$symbol %in% co_all_HALLMARK,]

co_nature09247_phf = pvalue_hr[pvalue_hr$symbol %in% co_nature09247,]

#################################################################  HALLMARK_INTERFERON_GAMMA_RESPONSE

# merge(x = df1, y = df2, by = "CustomerId", all.y=TRUE)

colnames(HALLMARK_INTERFERON_GAMMA_RESPONSE)[4] = "symbol"
co_all_HALLMARK_phf$gene_names = rownames(co_all_HALLMARK_phf)
co_all_HALLMARK_phf = co_all_HALLMARK_phf[,c("symbol","gene_names","pvalue_data","hr_data","fdr_data")]

HALLMARK_new = merge(x = HALLMARK_INTERFERON_GAMMA_RESPONSE, y = co_all_HALLMARK_phf, by = "symbol", all.x=TRUE)

# write.csv(HALLMARK_new,"D:/R_workspace/survival/PPT/methy_data/60data/outp/HALLMARK_new.csv")
# write.xlsx(HALLMARK_new,file= "D:/R_workspace/survival/PPT/methy_data/60data/outp/HALLMARK_new.xlsx",colNames = TRUE)
co_all_HALLMARK_feature = intersect(HALLMARK_INTERFERON_GAMMA_RESPONSE_GENE.SYMBOL,feature$names)


##################################################################  nature09247
colnames(co_nature09247_phf)[4] = "Symbol"
co_nature09247_phf$gene_names = rownames(co_nature09247_phf)
co_nature09247_phf = co_nature09247_phf[,c("Symbol","gene_names","pvalue_data","hr_data","fdr_data")]

nature09247_new = merge(x = nature09247, y = co_nature09247_phf, by = "Symbol", all.x=TRUE)

# write.csv(nature09247_new,"D:/R_workspace/survival/PPT/methy_data/60data/find_gene_data/nature09247_new.csv")

co_all_nature09247_feature = intersect(nature09247_GENE.SYMBOL,feature$names)


################################################


pvalue_hr = read.csv('D:/R_workspace/survival/PPT/methy_data/60data/combat_methy_pvalue_hr_fdr.csv',row.names = 1,stringsAsFactors = FALSE)

genes_data2 = read.csv('D:/R_workspace/survival/PPT/methy_data/multianno.txt',stringsAsFactors = FALSE,sep="")
gene_t1 = genes_data2[c(1,2,3,6,7)]
gene_t2 = unite(gene_t1, "index_methy", c(1,2,3))
gene_t2$count_nchar = apply(gene_t2["Gene.ensGene"],1,nchar)
gene_t3 = subset(gene_t2,count_nchar == 15)
gene_unique = unique(gene_t3$Gene.ensGene)

###############################################################   HALLMARK_INTERFERON_GAMMA_RESPONSE
names_HALLMARK_new = HALLMARK_new$gene_names
asa = str_split_fixed(names_HALLMARK_new,"_",2)


names_HALLMARK_new_data_frame = data.frame(asa)
names_HALLMARK_new_data_frame$gene_fullnames = as.character(HALLMARK_new$gene_names)
colnames(names_HALLMARK_new_data_frame)[1:2] = c("Gene.ensGene","af_genes")
names_HALLMARK_new_excu = as.character(names_HALLMARK_new_data_frame$X1)




co_gene = intersect(names_HALLMARK_new_excu,gene_unique)

gene_t4 = gene_t3[gene_t3$Gene.ensGene %in% co_gene,]


HALLMARK_new_me = merge(x = gene_t4, y = names_HALLMARK_new_data_frame, by = "Gene.ensGene", all.x=TRUE)


HALLMARK_new_me_t5 = HALLMARK_new_me[c(1,2,3,5,6)]

pvalue_hr_me = pvalue_hr

pvalue_hr_me$index_methy = rownames(pvalue_hr)

HALLMARK_new_me_phf = merge(x = HALLMARK_new_me_t5, y = pvalue_hr_me, by = "index_methy", all.x=TRUE)



index_al = c("index_methy","Gene.ensGene","af_genes","gene_fullnames","Func.ensGene","pvalue_data","hr_data","fdr_data")
HALLMARK_new_me_phf = HALLMARK_new_me_phf[,index_al]

# write.csv(HALLMARK_new_me_phf,"D:/R_workspace/survival/PPT/methy_data/60data/find_gene_data/HALLMARK_methy_phf.csv")


###############################################################   nature09247

pvalue_hr = read.csv('D:/R_workspace/survival/PPT/methy_data/60data/combat_methy_pvalue_hr_fdr.csv',row.names = 1,stringsAsFactors = FALSE)

genes_data2 = read.csv('D:/R_workspace/survival/PPT/methy_data/multianno.txt',stringsAsFactors = FALSE,sep="")
gene_t1 = genes_data2[c(1,2,3,6,7)]
gene_t2 = unite(gene_t1, "index_methy", c(1,2,3))
gene_t2$count_nchar = apply(gene_t2["Gene.ensGene"],1,nchar)
gene_t3 = subset(gene_t2,count_nchar == 15)
gene_unique = unique(gene_t3$Gene.ensGene)


names_nature09247_new = nature09247_new$gene_names
asa = str_split_fixed(names_nature09247_new,"_",2)


names_nature09247_new_data_frame = data.frame(asa)
names_nature09247_new_data_frame$gene_fullnames = as.character(nature09247_new$gene_names)
colnames(names_nature09247_new_data_frame)[1:2] = c("Gene.ensGene","af_genes")
names_nature09247_new_excu = as.character(names_nature09247_new_data_frame$Gene.ensGene)



co_gene = intersect(names_nature09247_new_excu,gene_unique)

gene_t4 = gene_t3[gene_t3$Gene.ensGene %in% co_gene,]


nature09247_new_me = merge(x = gene_t4, y = names_HALLMARK_new_data_frame, by = "Gene.ensGene", all.x=TRUE)


HALLMARK_new_me_t5 = HALLMARK_new_me[c(1,2,3,5,6)]

pvalue_hr_me = pvalue_hr

pvalue_hr_me$index_methy = rownames(pvalue_hr)

HALLMARK_new_me_phf = merge(x = HALLMARK_new_me_t5, y = pvalue_hr_me, by = "index_methy", all.x=TRUE)



index_al = c("index_methy","Gene.ensGene","af_genes","gene_fullnames","Func.ensGene","pvalue_data","hr_data","fdr_data")
HALLMARK_new_me_phf = HALLMARK_new_me_phf[,index_al]

# write.csv(HALLMARK_new_me_phf,"D:/R_workspace/survival/PPT/methy_data/60data/find_gene_data/HALLMARK_methy_phf.csv")


