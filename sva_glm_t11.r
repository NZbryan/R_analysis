require(sva)
require(ggplot2)
require(glmnet)
require(survival)
#require(reshape2)
require(rms)
require(caret)
#require(affy)
#require(affydata)

t1 = read.csv('survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/R_workspace/survival/expr_data.csv',row.names = 1)
rownames(t1) = colnames(t2)

#outlier_samples =c(A7,GX009,X2C2,X2D1)
t1 = t1[!rownames(t1) %in% c("A7","GX009","X2C2","X2D1","A8"),]
t2 = subset(t2,select = -c(A7,GX009,X2C2,X2D1,A8))

## rm all zero by the feature
rm_zero_func<-function(df,...)
{
  return(all(df==0))
}
col_t2 = apply(t2,1,rm_zero_func)
t2 = t2[!col_t2,]
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
#write.csv(combat_data_protocol,'D:/R_workspace/survival/PPT/combat_data_protocol.csv')
#write.csv(t1,'D:/R_workspace/survival/PPT/patient_info.csv')


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

#### strata
require(sampling)
type_0 = round(2/3*sum(t1$cencor.status==0))
type_1 = round(2/3*sum(t1$cencor.status==1))

#pdf('D:/R_workspace/survival/PPT/test1.pdf',width=15,height = 12)
#for(k in 1:10){
#num_set = as.integer(paste(20171120,k,sep=''))
#set.seed(num_set)
#set.seed(201711200)
sub = strata(t1,stratanames="cencor.status",size=c(type_1,type_0),method = "srswor")
bef_tr_data = scale_median_data[,sub$ID_unit]
bef_te_data = scale_median_data[,-sub$ID_unit]

adj_data = t1[sub$ID_unit,]
y_te_data = t1[-sub$ID_unit,]

mod = model.matrix(~as.factor(cencor.status), data=adj_data)
mod0 = model.matrix(~1, data = adj_data)

pValuesComBat = f.pvalue(bef_tr_data,mod,mod0)

qValuesComBat = p.adjust(pValuesComBat,method="BH")

pValuesComBat_005 = pValuesComBat[pValuesComBat<0.05]
pValuesComBat_005 = rownames(data.frame(pValuesComBat_005))

sample_data = as.data.frame(t(bef_tr_data))

pv = sample_data[pValuesComBat_005]
pv$gender = adj_data$gender

x = as.matrix(pv)
y = adj_data$survival.month



# cvfit_1 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=1,nlambda = 150)
# cvfit_2 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0,nlambda = 150)
# cvfit_3 = cv.glmnet(x, y, family = "gaussian",maxit=15000,alpha=0.2,nlambda = 150)
# par(mfrow=c(3,1))
# plot(cvfit_1)
# plot(cvfit_2)
# cvfit_3 = cv.glmnet(x, y, family = "gaussian",alpha=0.1,nlambda = 150)
# plot(cvfit_3)

# coefficients<-coef(cvfit_3,s=0.001)
# #coefficients[which(coefficients != 0)] 
# chosen_index = rownames(coefficients)[which(coefficients != 0)]
# 

fit_04 = glmnet(x, y, family="gaussian", alpha=0.42,lambda = 1.3750)

a11 = fit_04$beta


beta = a11[which(a11 != 0)]
chosen_names = rownames(a11)[which(a11 != 0)]
chosen_features = data.frame(chosen_names,beta)
#write.csv(chosen_features,'D:/R_workspace/survival/PPT/chosen_features_1.csv')



data_chosen_m = as.data.frame(x)
data_chosen = data_chosen_m[chosen_names]
data_chosen_x = as.matrix(data_chosen)
#fit_05 = glmnet(data_chosen_x, y, family="gaussian", alpha=0.42,lambda = 1.3750)

#### train & test set
bef_te_data = rbind(bef_te_data,y_te_data$gender)
rownames(bef_te_data)[dim(bef_te_data)[1]] = 'gender'

trainData_x = data_chosen_x
testData_x = t(bef_te_data[chosen_names,])
trainData_y = y
testData_y = y_te_data$survival.month

## predict
fit_a = glmnet(trainData_x, trainData_y, family="gaussian", alpha=0.42,lambda = 1.3750)
pre_y_a = predict(fit_a,testData_x)

data_a = data.frame(testData_y,pre_y_a)
colnames(data_a)=c('y_a','predict_a')

y_a_status = y_te_data$cencor.status
data_a$y_a_status = y_a_status

# log rank & km curve
data_a$predict_data = data_a$predict_a
md_va = data_a$predict_a>39
data_a$predict_a[md_va] = 1
data_a$predict_a[!md_va] = 2


SurvObj_a <- with(data_a, Surv(y_a, y_a_status == 1))
km_a = npsurv(SurvObj_a ~ predict_a, data = data_a)

log_rank_predict_a=survdiff(SurvObj_a ~ predict_a, data = data_a)
log_rank_predict_a_pvalue = pchisq(log_rank_predict_a$chisq, length(log_rank_predict_a$n)-1, lower.tail = FALSE)

survplot(km_a,xlab = 'months', n.risk   = TRUE,cex.n.risk = 1,type="kaplan-meier",lty=1,lwd=2,
         label.curves = list(method = "arrow", cex = 1.2),col=c("orange","blue"),abbrev.label=TRUE)
text(60,1,sprintf("p-value of log-rank test: %.10f", log_rank_predict_a_pvalue),col='brown2')
title("high-risk:2   low-risk:1")




