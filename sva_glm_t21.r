require(sva)
require(ggplot2)
require(glmnet)
require(survival)
#require(reshape2)
require(rms)
require(caret)
#require(affy)
#require(affydata)

t1 = read.csv('D:/python_workspace/survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/python_workspace/survival/expression.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')

rownames(t1) = colnames(t2)
#outlier_samples =c(A7,GX009,X2C2,X2D1)
t1 = t1[!rownames(t1) %in% c("A7","GX009","X2C2","X2D1","A8"),]
t2 = subset(t2,select = -c(A7,GX009,X2C2,X2D1,A8))
patient_info = t1
patient_info$SurvObj <- with(patient_info, Surv(survival.month, cencor.status == 1))

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

t3 = log2(t2)
t3[t3==-Inf]=0
col_t3 = apply(t3,1,rm_zero_func)
t3 = t3[!col_t3,]

combat_expr_data  = ComBat(dat=t3, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)
combat_expr_data = t(combat_expr_data)
colnames_expr = rownames(combat_expr_data)

pvalue_length = dim(combat_expr_data)[2]
pvalue_data = array(1:pvalue_length)
hr_data = array(1:pvalue_length)
for (i in 1:pvalue_length)
{
  model = coxph(SurvObj~combat_expr_data[,i],data = patient_info)
  pra = glance(model)$p.value.log
  pvalue_data[i] = pra
  hr_data[i] = model$coefficients
}

# pvalue_logrank=as.data.frame(pvalue_data)
# hr_logrank = as.data.frame(hr_data)
hr_data = exp(hr_data)
pvalue_hr = data.frame(pvalue_data,hr_data)
rownames(pvalue_hr) = colnames_expr


b11=sort(pvalue_hr$hr_data,decreasing = TRUE)[1000]
pvalue_hr_005_rank1000 = subset(pvalue_hr,hr_data>b11 & pvalue_data<0.05)
pvalue_hr_005_rank1000_names = rownames(pvalue_hr_005_rank1000)

x = subset(combat_expr_data,select = pvalue_hr_005_rank1000_names)
x = as.data.frame(x)
x$gender = t1$gender
x = as.matrix(x)
y = patient_info$SurvObj


set.seed(2017118009)
cvfit_cox = cv.glmnet(x, y, family = "cox",alpha=0.5,nlambda = 150)
plot(cvfit_cox)

fit_cox = glmnet(x, y, family="cox",alpha=0.5,lambda = cvfit_cox$lambda.min)


## predict

# y_predict = predict(fit_cox,x,type = "response")
y_predict = predict(fit_cox,x)

data_a = data.frame(y,y_predict)
colnames(data_a)=c('SurvObj','y_predict')

# log rank & km curve
data_a$signature = data_a$y_predict
md_va = data_a$signature>median(y_predict)
data_a$signature[md_va] = "High-risk"
data_a$signature[!md_va] = "Low-risk"
# status_mid1 = table(data_a$predict_binary,data_a$y_true_status)[1]+
#   table(data_a$predict_binary,data_a$y_true_status)[4]
# status_mid2 = table(data_a$predict_binary,data_a$y_true_status)[2]+
#   table(data_a$predict_binary,data_a$y_true_status)[3]
# if(status_mid1<status_mid2)
# {
#   data_a$predict_binary[data_a$predict_binary == 1] = 2
#   data_a$predict_binary[data_a$predict_binary == 0] = 1
#   data_a$predict_binary[data_a$predict_binary == 2] = 0
# }
# 
if(length(unique(data_a$signature))==1){next}

km_a = npsurv(SurvObj ~ signature, data = data_a)

log_rank_predict_a=survdiff(SurvObj ~ signature, data = data_a)
log_rank_predict_a_pvalue = pchisq(log_rank_predict_a$chisq, length(log_rank_predict_a$n)-1, lower.tail = FALSE)

par(cex=1)
survplot(km_a,xlab = 'months', n.risk   = TRUE,y.n.risk=0,cex.n.risk = 1,type="kaplan-meier",lty=1,lwd=3,
         label.curves = list(method = "arrow", cex = 1.2),col=c("orange","blue"),abbrev.label=FALSE)
text(60,0.65,sprintf("p-value of log-rank test: %.10f", log_rank_predict_a_pvalue),col='brown2')
title_K = paste("test:",k)
title(title_K)






