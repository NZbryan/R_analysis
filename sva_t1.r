require(sva)
require(pamr)

t1 = read.csv('D:/python_workspace/survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/python_workspace/survival/expression.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
#data_1 = t1[c('gender','age','survival.month','cencor.status')]
#t2 = as.matrix(t2)
#t3 = t(t2)
#t3 = as.data.frame(t3)
#rownames(t1) = colnames(t2)
#n.sv = num.sv(t2,mod,method="leek")
#svobj = sva(t2,mod,mod0,n.sv=n.sv)

t6 = read.csv('D:/R_workspace/survival/t6.csv',row.names = 1)

mod = model.matrix(~as.factor(gender), data=t1)
mod0 = model.matrix(~1,data=t1)



## sva
t6 = as.matrix(t6)
nsv6 = num.sv(t6,mod,method="leek")
svobj6 = sva(t6,mod,mod0,n.sv=10)

pValues = f.pvalue(t6,mod,mod0)
qValues = p.adjust(pValues,method="BH")

modSv = cbind(mod,svobj6$sv)
mod0Sv = cbind(mod0,svobj6$sv)
pValuesSv = f.pvalue(t6,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

## limma package
fit = lmFit(t6,modSv)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,9)),"C2"=c(0,-1,1,rep(0,9)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")

### ComBat
batch = t1$china2008.stage
modcombat = model.matrix(~1, data=t1)
combat_edata = ComBat(dat=t6, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

## fsva (remove latent variable)
set.seed(12354)
trainIndicator = sample(1:76,size=35,replace=F)
testIndicator = (1:76)[-trainIndicator]
trainData = t6[,trainIndicator]
testData = t6[,testIndicator]
trainPheno = t1[trainIndicator,]
testPheno = t1[testIndicator,]

mydata = list(x=trainData,y=trainPheno$cencor.status)
mytrain = pamr.train(mydata)

table(pamr.predict(mytrain,testData,threshold=2),testPheno$cencor.status)



trainMod = model.matrix(~cencor.status,data=trainPheno)
trainMod0 = model.matrix(~1,data=trainPheno)
trainSv = sva(trainData,trainMod,trainMod0)



count_0 = t(trainData)

f<-function(x)
{
  sum(x==0)
}

data_0 = apply(count_0,2,f)

index_11 = names(data_0[(data_0!=35)])

trainData = subset(trainData,select = index_11)









