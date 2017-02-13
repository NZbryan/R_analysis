library(glmnet)
library(rms)
library(survival)
t1 = read.csv('D:/python_workspace/survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/python_workspace/survival/expression.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
#data_1 = t1[c('gender','age','survival.month','cencor.status')]
t6 = read.csv('D:/R_workspace/survival/t6.csv',row.names = 1)
t6 = t(t6)
#t6 = as.data.frame(t6)

t3 = t(t2)
t3 = as.data.frame(t3)
#t3$sequencing.protocol = t1$sequencing.protocol


t1$gender[which(t1$gender=='M')]=1
t1$gender[which(t1$gender=='F')]=0

#t1$SurvObj <- with(t1, Surv(survival.month, cencor.status == 1))
surv<-Surv(t1$survival.month,t1$cencor.status)

t3$gender = t1$gender

x<-as.matrix(t3)

cv.fit<-cv.glmnet(x,surv,family="cox",maxit=15000,alpha=1,nlambda = 150)


##  choose features
coefficients<-coef(cv.fit,s=0.003)

coefficients[which(rownames(coefficients)=='gender')]
#coefficients[which(coefficients != 0)] 
chosen_index = rownames(coefficients)[which(coefficients != 0)]



t4 = t3[chosen_index]
x_new = as.matrix(t4)
new_fit<-cv.glmnet(x_new,surv,family="cox")


plot(fit, xvar="lambda", label=TRUE)


predict(new_fit, newx=x_new, type="response", s=new_fit$lambda.min)
fit = glmnet(x,surv,family="cox",maxit=15000,alpha=1)
predict(fit, newx=x, type="response", s=0.003)






