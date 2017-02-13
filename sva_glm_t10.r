require(sva)
require(ggplot2)
require(glmnet)
require(survival)
#require(reshape2)
require(rms)
require(caret)
#require(affy)
#require(affydata)
require(ggfortify)
require(data.table)

t1 = read.csv('survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/R_workspace/survival/expr_data.csv',row.names = 1)
rownames(t1) = colnames(t2)

#outlier_samples =c(A7,GX009,X2C2,X2D1)
# t1 = t1[!rownames(t1) %in% c("A7","GX009","X2C2","X2D1","A8"),]
# t2 = subset(t2,select = -c(A7,GX009,X2C2,X2D1,A8))

## rm all zero by the feature
rm_zero_func<-function(df,...)
{
  return(all(df==0))
}
col_t2 = apply(t2,1,rm_zero_func)
t2 = t2[!col_t2,]
t2 = as.matrix(t2)
t2 = t(t2)


t1$gender[which(t1$gender=='F')]=1
t1$gender[which(t1$gender=='M')]=2
t1$gender = as.integer(t1$gender)

# sequencing.protocol  batch
batch = t1$sequencing.protocol

batch[batch=='anchordx' & t1$gender==1] = "anchordx_F"
batch[batch=='anchordx' & t1$gender==2] = "anchordx_M"
batch[batch=='novogene' & t1$gender==1] = "novogene_F"
batch[batch=='novogene' & t1$gender==2] = "novogene_M"

cor_da = data.frame(batch)
names(cor_da) = "Species"



#png("D:/R_workspace/survival/PPT/nm1.png",width = 1706,height = 960)
autoplot(prcomp(t2), colour = "Species",data = cor_da,size=5,frame = TRUE,frame.type = 'norm')
text('normal data')

#autoplot(prcomp(t2), colour = c("red","green","orange","blue"))
#dev.off()

text(-0.2,0.4,sprintf("anchordx & F: red\nanchordx & M: green\nnovogene & F: orange\nnovogene & M: blue"))


t2 = t(t2)
modcombat = model.matrix(~1, data = t1)
combat_data_protocol  = ComBat(dat=t2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

combat_data_protocol = t(combat_data_protocol)
#png("D:/R_workspace/survival/PPT/nm1_combat.png",width = 1706,height = 960)

autoplot(prcomp(combat_data_protocol), colour = "Species",data = cor_da,size=5,frame = TRUE,frame.type = 'norm')
#dev.off()


## rm A7,GX009,X2C2,X2D1
#autoplot(prcomp(combat_data_protocol), colour = "Species",data = cor_da,label = TRUE)

## PCA

PCA_original = prcomp(t2)

eig <- (PCA_original$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)

plot_pca_data = eig.decathlon2.active[1:10,]
plot_pca_data_names = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

barplot(plot_pca_data[, 2], names.arg=plot_pca_data_names, ylim=c(0,60),
        main = "Principal Component Analysis (original data)",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")

lines(x = 1:nrow(plot_pca_data), 
      plot_pca_data[, 2], 
      type="b", pch=19, col = "red")

screeplot(PCA_original)



###

combat_expr_data = fread('D:/R_workspace/survival/PPT/combat_expr_data.csv')
patient_info = read.csv('D:/R_workspace/survival/PPT/patient_info.csv',stringsAsFactors = FALSE,row.names = 1)

combat_expr_data = as.data.frame(combat_expr_data)
rownames(combat_expr_data) = combat_expr_data$V1
combat_expr_data = combat_expr_data[,-c(1)]


patient_info$gender[which(patient_info$gender=='F')]=1
patient_info$gender[which(patient_info$gender=='M')]=2
patient_info$gender = as.integer(patient_info$gender)

# sequencing.protocol  batch
batch = patient_info$sequencing.protocol

batch[batch=='anchordx' & patient_info$gender==1] = "anchordx_F"
batch[batch=='anchordx' & patient_info$gender==2] = "anchordx_M"
batch[batch=='novogene' & patient_info$gender==1] = "novogene_F"
batch[batch=='novogene' & patient_info$gender==2] = "novogene_M"

cor_da = data.frame(batch)
names(cor_da) = "Species"
#PLOT pca
autoplot(prcomp(combat_expr_data), colour = "Species",data = cor_da,size=5,frame = TRUE,frame.type = 'norm')
#dev.off()
