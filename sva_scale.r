require(sva)
require(ggplot2)
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




#########################



col_batch = batch
col_batch[col_batch==1] = 'brown1'
col_batch[col_batch==2] = 'green'

#normal_data
bef_data = as.data.frame(t2)
bef_data[bef_data==0]=NA
boxplot(bef_data)


#combat_data
aft_data = combat_edata
aft_data[aft_data==0]=NA

#median——scale
j2 = log(aft_data)
mean_median = mean(apply(aft_data,2,median,na.rm = TRUE),na.rm = TRUE)

scale_median_2 <- function(df){
  a = median(df,na.rm = TRUE)
  b = sd(df,na.rm = TRUE)
  c = (df-a+mean_median)/b
  return(c)
}
scale_median_combat = apply(j2,2,scale_median_2)



# png("survival/data_expr//combat_scale.png",width = 2560,height = 1440)
par(mfrow = c(3,1))

par(cex=1.2)
boxplot(log(bef_data),col=c(col_batch),main="original_data")
boxplot(log(aft_data),col=c(col_batch),main="combat_data")
boxplot(scale_median_combat,col=c(col_batch),main="scale_median")
# dev.off()




scale_median_3 <- function(df){
  a = median(df,na.rm = TRUE)
  c = (df-a+mean_median)
  return(c)
}

scale_median_data =  apply(aft_data,2,scale_median_3)
a11 = apply(scale_median_data,2,median,na.rm = TRUE)





