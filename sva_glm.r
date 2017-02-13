require(sva)
require(ggplot2)
require(affy)
require(affydata)

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

## median scale
aft_data = combat_edata
aft_data[aft_data==0]=NA
mean_median = mean(apply(aft_data,2,median,na.rm = TRUE),na.rm = TRUE)

scale_median_3 <- function(df){
  a = median(df,na.rm = TRUE)
  c = (df-a+mean_median)
  return(c)
}

scale_median_data =  apply(combat_edata,2,scale_median_3)
#a11 = apply(scale_median_data,2,median,na.rm = TRUE)



mod = model.matrix(~as.factor(cencor.status), data=t1)
mod0 = model.matrix(~1, data = t1)

pValuesComBat = f.pvalue(scale_median_data,mod,mod0)

qValuesComBat = p.adjust(pValuesComBat,method="BH")

pValuesComBat_005 = pValuesComBat[pValuesComBat<0.05]
pValuesComBat_005 = rownames(data.frame(pValuesComBat_005))

sample_data = as.data.frame(t(scale_median_data))

pv = sample_data[pValuesComBat_005]
pv$gender = t1$gender

x = as.matrix(pv)
y = t1$survival.month




cvfit_1 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=1,nlambda = 150)

cvfit_2 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0,nlambda = 150)

cvfit_3 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0.5,nlambda = 150)


png("survival/data_expr/lasso_ridge_elastic.png",width = 1680,height = 1050)
par(mfrow=c(1,3))
plot(cvfit_1, main="LASSO")
plot(cvfit_2, main="Ridge")
plot(cvfit_3, main="Elastic Net")
dev.off()


cvfit = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0.45,nlambda = 150)
plot(cvfit)

fit_04 = glmnet(x, y, family="gaussian", alpha=0.5,lambda = 3.195319)



x_test = t(t2)


cvfit_1 = cv.glmnet(x_test, y, family = "gaussian",maxit=15000,alpha=1,nlambda = 150)
cvfit_2 = cv.glmnet(x_test, y, family = "gaussian",maxit=15000,alpha=0,nlambda = 150)
cvfit_3 = cv.glmnet(x_test, y, family = "gaussian",maxit=15000,alpha=0.5,nlambda = 150)
png("survival/data_expr/lasso_ridge_elastic_bef.png",width = 1680,height = 1050)
par(mfrow=c(1,3))
plot(cvfit_1, main="LASSO")
plot(cvfit_2, main="Ridge")
plot(cvfit_3, main="Elastic Net")
dev.off()