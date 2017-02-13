require(sva)
require(ggplot2)
#require(affy)
#require(affydata)

t1 = read.csv('survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/R_workspace/survival/expr_data.csv',row.names = 1)
t2 = as.matrix(t2)

t1$gender[which(t1$gender=='F')]=1
t1$gender[which(t1$gender=='M')]=2
t1$gender = as.integer(t1$gender)

# sequencing.protocol  batch
batch = t1$sequencing.protocol

batch[batch=='anchordx' & t1$gender==1] = 1
batch[batch=='anchordx' & t1$gender==2] = 2
batch[batch=='novogene' & t1$gender==1] = 3
batch[batch=='novogene' & t1$gender==2] = 4

modcombat = model.matrix(~1, data = t1)

combat_data_protocol  = ComBat(dat=t2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

aft_data = combat_data_protocol
aft_data[aft_data==0]=NA
mean_median = mean(apply(aft_data,2,median,na.rm = TRUE),na.rm = TRUE)

scale_median_3 <- function(df){
  a = median(df,na.rm = TRUE)
  c = (df-a+mean_median)
  return(c)
}

scale_median_data =  apply(combat_data_protocol,2,scale_median_3)
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



set.seed(123589)
cvfit_1 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=1,nlambda = 150)

cvfit_2 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0,nlambda = 150)

set.seed(12356896)
cvfit_3 = cv.glmnet(x, y, family = "poisson",maxit=15000,alpha=0.42,nlambda = 150)
plot(cvfit_3)

set.seed(12356896)
cvfit_4 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0.42,nlambda = 150)
plot(cvfit_4)

coefficients<-coef(cvfit_3,s=1.5)
#coefficients[which(coefficients != 0)] 
chosen_index = rownames(coefficients)[which(coefficients != 0)]

fit_04 = glmnet(x, y, family="gaussian", alpha=0.42,lambda = 1.3750)

a11 = fit_04$beta


beta = a11[which(a11 != 0)]
chosen_names = rownames(a11)[which(a11 != 0)]
chosen_features = data.frame(chosen_names,beta)
#write.csv(chosen_features,'D:/R_workspace/survival/PPT/chosen_features_1.csv')

pre = predict(fit_04,x)


data_1 = data.frame(y,pre)
colnames(data_1)=c('y','predict')

p = ggplot(data = data_1, aes(x = y, y = predict))+ geom_point(col='blue')
p + geom_abline(intercept = 0, slope = 1,col = 'red',size=1.5)+geom_smooth(method=lm)




par(mfrow=c(2,1))
ggplot(data = data_1, aes(x = y))+geom_histogram(binwidth=3,color='orange')+geom_vline(aes(xintercept=39), col='blue', linetype="dashed",size=2)
ggplot(data = data_1, aes(x = predict))+geom_histogram(binwidth=3,color='blue')+geom_vline(aes(xintercept=39), col='blue', linetype="dashed",size=2)










