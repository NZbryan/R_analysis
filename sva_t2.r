require(sva)
require(ggplot2)
require(affy)
require(affydata)

t1 = read.csv('survival/clinical.info.data.frame.txt',header = T,row.names = 1,stringsAsFactors = FALSE,sep='')
t2 = read.csv('D:/R_workspace/survival/t6.csv',row.names = 1)
batch
t2 = as.matrix(t2)
hist(log2(t2),breaks=100,col="blue") 

batch = t1$sequencing.protocol

batch[batch=='anchordx'] = 1
batch[batch=='novogene'] = 2

modcombat = model.matrix(~1, data = t1)

combat_edata = ComBat(dat=t2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)


par(mfrow=c(1,2))
hist(log2(expr_data),breaks=100,col="blue") 
hist(log2(combat_edata),breaks=100,col="blue") 


mod = model.matrix(~as.factor(cencor.status), data=t1)
mod0 = model.matrix(~1, data = t1)

pValuesComBat = f.pvalue(combat_edata,mod,mod0)

qValuesComBat = p.adjust(pValuesComBat,method="BH")



length(pValuesComBat[pValuesComBat<0.05])/length(pValuesComBat)
pValues = f.pvalue(t2,mod,mod0)
length(pValues[pValues<0.05])/length(pValues)



pValuesComBat = data.frame(pValuesComBat)
pValues = data.frame(pValues)
pValuesComBat$type = 'pValuesComBat'
pValues$type = 'pValues'

names(pValues)[1] = names(pValuesComBat)[1]
co_pvalue = rbind(pValuesComBat,pValues)
names(co_pvalue)[1] = 'values'
ggplot(co_pvalue)+aes(x=values)+geom_histogram(position = 'dodge',binwidth=0.01)+aes(fill=type)

expr_data = t2
par(mfrow=c(1,2))
hist(log2(expr_data),breaks=100,col="blue") 
hist(log2(combat_edata),breaks=100,col="blue") 



affybatch <- ReadAffy(celfile.path = "survival/data_expr/GSE7890_RAW.tar" )

AffyData<-ReadAffy(widget=TRUE) 
exprsSet = rma(AffyData)


write.csv(combat_edata,file='survival/combat_edata')


data(Dilution)
eset <- rma(Dilution)

exp.RMA<-exprs(eset)

par(mfrow=c(1,2))
boxplot(Dilution,col="red")
boxplot(data.frame(exp.RMA),col="blue")



