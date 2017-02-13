require(sva)
require(ggplot2)
require(glmnet)
# require(affy)
# require(affydata)

t1 = read.csv('survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/R_workspace/survival/expr_data.csv',row.names = 1)
t2 = as.matrix(t2)

t1$gender[which(t1$gender=='F')]=1
t1$gender[which(t1$gender=='M')]=2
t1$gender = as.integer(t1$gender)

# combat
batch = t1$sequencing.protocol

batch[batch=='anchordx'] = 1
batch[batch=='novogene'] = 2

modcombat = model.matrix(~1, data = t1)

combat_edata = ComBat(dat=t2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

mod = model.matrix(~as.factor(cencor.status), data=t1)
mod0 = model.matrix(~1, data = t1)

pValuesComBat = f.pvalue(combat_edata,mod,mod0)

qValuesComBat = p.adjust(pValuesComBat,method="BH")

pValuesComBat_005 = pValuesComBat[pValuesComBat<0.05]
pValuesComBat_005 = rownames(data.frame(pValuesComBat_005))

sample_data = as.data.frame(t(combat_edata))

pv = sample_data[pValuesComBat_005]
pv$gender = t1$gender

x = as.matrix(pv)
y = t1$survival.month





cvfit_1 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=1,nlambda = 150)

cvfit_2 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0,nlambda = 150)

cvfit_3 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0.3,nlambda = 150)
plot(cvfit_3)
# png("survival/lasso_ridge_elastic.png",width = 1680,height = 1050)
par(mfrow=c(1,3))
plot(cvfit_1, main="LASSO")
plot(cvfit_2, main="Ridge")
plot(cvfit_3, main="Elastic Net")
# dev.off()

par(mar=c(5,4,7,2)+0.1,cex=1.2)
plot(cvfit_3)
title("Elastic-Net")

fit = glmnet(x, y, family="gaussian", alpha=0.8,lambda = 1.111611)

pre =  predict(

fit_11 = glm()


cvfit = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0.4,nlambda = 150)

fit_04 = glmnet(x, y, family="gaussian", alpha=0.4,lambda = 2.846856)



data_1 = pv


data_1$y = y
data_1 = as.data.frame(data_1)

b11 = glm(y~., family = gaussian, data = data_1)


p<-ggplot(data=t2, aes(x=t2,y=Count))+geom_boxplot(aes(fill=Organ_type))



bef_data = t(t2)
bef_data = as.data.frame(bef_data)
batch = as.integer(batch)
bef_data$type = batch

p<-ggplot(data=Data, aes(x=Expression_level,y=Count))+geom_boxplot(aes(fill=Organ_type))

require(ggplot2)

bef_data = as.data.frame(t2)
p<-ggplot(data=bef_data, aes(x=bef_data))+geom_boxplot()


bef_data[bef_data==0]=NA
boxplot(bef_data)
boxplot(log(bef_data),col='red')

aft_data = combat_edata
aft_data[aft_data==0]=NA
boxplot(log(bef_data),col='green')


scale_median_1 <- function(df){
  a = median(df,na.rm = TRUE)
  b = sd(df,na.rm = TRUE)
  c = (df-a)/b
  return(c)
}



j1 = combat_edata
j1[j1==0] = NA

j2 = log(j1)
scale_median_combat = apply(j2,2,scale_median_1)


boxplot(scale_median_combat,col='green')


scale_median_2 <- function(df){
  a = median(df,na.rm = TRUE)
  b = sd(df,na.rm = TRUE)
  c = (df-a+1)/b
  return(c)
}


j1 = combat_edata
j1[j1==0] = NA


j3 = apply(j2,2,median)


j2 = log(j1)
scale_median_combat = apply(j2,2,scale_median_2)


scale_median_combat = apply(j2,2,scale_median_2)
png("survival/combat_scale.png",width = 2560,height = 1440)
par(mfrow = c(3,1))
boxplot(log(bef_data),col='green',main="normal_data")
boxplot(log(j1),col='red',main="combat_data")
boxplot(scale_median_combat,col='orange',main="scale_median")
dev.off()


bef_data = as.data.frame(t2)
bef_data[bef_data==0]=NA

j1 = combat_edata
j1[j1==0] = NA

