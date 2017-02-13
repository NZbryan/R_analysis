require(sva)
require(data.table)
require(ggplot2)
require(glmnet)
require(survival)
#require(reshape2)
require(rms)
require(caret)
require(broom)
#require(affy)
#require(affydata)
options(warn=-1) 

combat_expr_data = fread('D:/R_workspace/survival/PPT/combat_expr_data.csv')
patient_info = read.csv('D:/R_workspace/survival/PPT/patient_info.csv',row.names = 1)
pvalue_logrank = read.csv('D:/R_workspace/survival/PPT/pvalue_logrank.csv',row.names = 1)

combat_expr_data = as.data.frame(combat_expr_data)
rownames(combat_expr_data) = combat_expr_data$V1
combat_expr_data = combat_expr_data[,-c(1)]
colnames_expr = colnames(combat_expr_data)

patient_info$SurvObj <- with(patient_info, Surv(survival.month, cencor.status == 1))

# pvalue_data = pvalue_logrank$pvalue_data
# padjust_logrank = p.adjust(pvalue_data,method="fdr")
# padjust_logrank = as.data.frame(padjust_logrank)
# rownames(padjust_logrank) = rownames(pvalue_logrank)

pvalue_logrank_005 = subset(pvalue_logrank,pvalue_data<0.05)
pvalue_logrank_005_names = rownames(pvalue_logrank_005)

x = subset(combat_expr_data,select = pvalue_logrank_005_names)
x$gender = patient_info$gender
x = as.matrix(x)
y = patient_info$survival.month


# set.seed(201711202)
# fit_cox = glmnet(x, patient_info$SurvObj, family="cox", alpha=0.1,lambda = 0.4)
# pred = predict(fit_cox,x)
# set.seed(201711600)
# cv_Folds = createFolds(1:dim(x)[1],k=5)
# #png("D:/R_workspace/survival/PPT/cv5.png",width = 2560,height = 1440)
# par(mfrow = c(5,1))

#num_set = as.integer(paste(20171120,k,sep=''))
#set.seed(num_set)
set.seed(201711202)
require(sampling)
type_0 = round(1/2*sum(patient_info$cencor.status==0))
type_1 = round(1/2*sum(patient_info$cencor.status==1))

sub = strata(patient_info,stratanames="cencor.status",size=c(type_1,type_0),method = "srswor")

test_x = x[sub$ID_unit,]
train_x = x[-sub$ID_unit,]

test_y = patient_info$SurvObj[sub$ID_unit,]
train_y = patient_info$SurvObj[-sub$ID_unit,]


## for
# for(k in 1:length(cv_Folds))
# {
#   
#   test_y = patient_info[cv_Folds[[k]],]$SurvObj
#   train_y = patient_info[-cv_Folds[[k]],]$SurvObj
#   
#   test_x = x[cv_Folds[[k]],]
#   train_x = x[-cv_Folds[[k]],]

cvfit_cox = cv.glmnet(x, patient_info$SurvObj, family = "cox",alpha=0.1,nlambda = 150)

## train the model
fit_cox = glmnet(train_x, train_y, family="cox",alpha=0.1,lambda = cvfit_cox$lambda.min)

## predict

y_predict = predict(fit_cox,test_x,type = "response")


data_a = data.frame(test_y,y_predict)
colnames(data_a)=c('SurvObj','y_predict')

# log rank & km curve
data_a$predict_binary = data_a$y_predict
md_va = data_a$predict_binary>2.5
data_a$predict_binary[md_va] = 1
data_a$predict_binary[!md_va] = 0
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
km_a = npsurv(SurvObj ~ predict_binary, data = data_a)

log_rank_predict_a=survdiff(SurvObj ~ predict_binary, data = data_a)
log_rank_predict_a_pvalue = pchisq(log_rank_predict_a$chisq, length(log_rank_predict_a$n)-1, lower.tail = FALSE)


#png("D:/R_workspace/survival/PPT/strata_1.png",width = 1706,height = 960)
par(mfrow = c(2,1))
survplot(km_a,xlab = 'months', n.risk   = TRUE,cex.n.risk = 1,type="kaplan-meier",lty=1,lwd=2,
       label.curves = list(method = "arrow", cex = 1.2),col=c("orange","blue"),abbrev.label=TRUE)
text(60,1,sprintf("p-value of log-rank test: %.10f", log_rank_predict_a_pvalue),col='brown2')
title_K = paste("test")
title(title_K)



#### 

test_x = x[-sub$ID_unit,]
train_x = x[sub$ID_unit,]

test_y = patient_info$SurvObj[-sub$ID_unit,]
train_y = patient_info$SurvObj[sub$ID_unit,]

cvfit_cox = cv.glmnet(train_x, train_y, family = "cox",alpha=0.1,nlambda = 150)

fit_cox = glmnet(train_x, train_y,family="cox",alpha=0.1,lambda = cvfit_cox$lambda.min)

## predict

y_predict = predict(fit_cox,test_x,type="response")


data_a = data.frame(test_y,y_predict)
colnames(data_a)=c('SurvObj','y_predict')

# log rank & km curve
data_a$predict_binary = data_a$y_predict
md_va = data_a$predict_binary>2.5
data_a$predict_binary[md_va] = 1
data_a$predict_binary[!md_va] = 0
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
km_a = npsurv(SurvObj ~ predict_binary, data = data_a)

log_rank_predict_a=survdiff(SurvObj ~ predict_binary, data = data_a)
log_rank_predict_a_pvalue = pchisq(log_rank_predict_a$chisq, length(log_rank_predict_a$n)-1, lower.tail = FALSE)

survplot(km_a,xlab = 'months', n.risk   = TRUE,cex.n.risk = 1,type="kaplan-meier",lty=1,lwd=2,
         label.curves = list(method = "arrow", cex = 1.2),col=c("orange","blue"),abbrev.label=TRUE)
text(60,1,sprintf("p-value of log-rank test: %.10f", log_rank_predict_a_pvalue),col='brown2')
title_K = paste("test")
title(title_K)
#dev.off()
