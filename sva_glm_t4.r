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





set.seed(12356896)
cvfit_4 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0.42,nlambda = 150)
plot(cvfit_4)

#chosen_index = rownames(coefficients)[which(coefficients != 0)]

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



data_chosen_m = as.data.frame(x)
data_chosen = data_chosen_m[chosen_names]

data_chosen_x = as.matrix(data_chosen)
fit_05 = glmnet(data_chosen_x, y, family="gaussian", alpha=0.42,lambda = 1.3750)




trainIndicator = sample(1:76,size=38,replace=F)

testIndicator = (1:76)[-trainIndicator]

trainData_x = data_chosen_x[trainIndicator,]

testData_x = data_chosen_x[testIndicator,]

trainData_y = y[trainIndicator]

testData_y = y[testIndicator]

fit_06 = glmnet(trainData_x, trainData_y, family="gaussian", alpha=0.42,lambda = 1.3750)

pre_y_06 = predict(fit_06,testData_x)




data_2 = data.frame(testData_y,pre_y_06)
colnames(data_2)=c('y','predict')

p = ggplot(data = data_2, aes(x = y, y = predict))+ geom_point(col='blue')
p + geom_abline(intercept = 0, slope = 1,col = 'red',size=1.5)+geom_smooth(method=lm)

cvfit_4 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0.42,nlambda = 150)

fit_07 = cv.glmnet(data_chosen_x, y, family = "gaussian",alpha=0.42)



fit_08 = glmnet(trainData_x, trainData_y, family="poisson", alpha=0.42,lambda = 1.3750)
pre_y_08 = predict(fit_08,testData_x)

data_3 = data.frame(testData_y,pre_y_08)
colnames(data_3)=c('y','predict')

p = ggplot(data = data_3, aes(x = y, y = predict))+ geom_point(col='blue')
p + geom_abline(intercept = 0, slope = 1,col = 'red',size=1.5)+geom_smooth(method=lm)


