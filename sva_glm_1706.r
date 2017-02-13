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

feature = read.xlsx('D:/R_workspace/survival/PPT/methy_data/RNASep.xlsx')


hclust_name = function(df)
{ 
  s1 = str_split_fixed(df,"_",2)[1]
}

b2 = apply(feature,1,hclust_name)

feature$names = b2

# unique(genes_data2$Gene.ensGene)
gene_t1 = genes_data2[c(1,2,3,7)]
gene_t2 = unite(gene_t1, "index_methy", c(1,2,3))
gene_t2$count_nchar = apply(gene_t2["Gene.ensGene"],1,nchar)
gene_t3 = subset(gene_t2,count_nchar == 15)


gene_unique = unique(gene_t3$Gene.ensGene)
co_gene = intersect(b2,gene_unique)

gene_t4 = gene_t3[gene_t3$Gene.ensGene %in% co_gene,]

gene_t5 = gene_t4[c(1,2)]

#write.csv(feature,"D:/R_workspace/survival/PPT/methy_data/60data/feature.csv")

pvalue_hr = read.csv('D:/R_workspace/survival/PPT/pvalue_hr_new.csv',row.names = 1)
wnt_genes = read.csv('D:/R_workspace/survival/PPT/methy_data/wnt_genes.csv', header = FALSE,row.names = 1,stringsAsFactors = FALSE)


y1 = rownames(pvalue_hr)
y2 = wnt_genes$V2



room1 = array(NA,length(y2))
find_genes = array(NA,length(y2))
for(i in 1:length(y2))
{
  mid1 = grep(paste0("_",y2[i],"$"),y1)
  if(length(mid1)>0)
  {
    room1[i]=y2[i]
    find_genes[i] = y1[mid1]
  }
  
}


fdr_data = p.adjust(pvalue_hr$pvalue_data, method="fdr")

pvalue_hr$fdr_data = fdr_data

pvalue_hr_fdr = pvalue_hr



find_genes_pvalue_hr_fdr = pvalue_hr_fdr[find_genes,]

find_genes_pvalue_hr_fdr$wnt_signaling_pathway = y2

#write.csv(find_genes_pvalue_hr_fdr,"D:/R_workspace/survival/PPT/methy_data/60data/find_genes_pvalue_hr_fdr.csv")

