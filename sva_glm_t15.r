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

pvalue_logrank_005 = subset(pvalue_logrank,pvalue_data<0.005)
pvalue_logrank_005_names = rownames(pvalue_logrank_005)

x = subset(combat_expr_data,select = pvalue_logrank_005_names)
x$gender = patient_info$gender
x = as.matrix(x)
y = patient_info$SurvObj



cvfit_cox = cv.glmnet(x, y, family = "cox",alpha=0.1,nlambda = 150)
plot(cvfit_cox)


# set.seed(201711600)
for(i in 0:9)
{
  sesd = as.integer(paste(201711800,i,sep = ""))
  set.seed(sesd)
  set.seed(2017118009)
  cv_Folds = createFolds(1:dim(x)[1],k=3)
  
  png_f = paste("D:/R_workspace/survival/PPT/cv_1.5/",i,".png",sep = "")
  #png("D:/R_workspace/survival/PPT/cv_Folds3_2.png",width = 650,height = 850)
  png(png_f,width = 1200,height = 1400)
  par(mfrow = c(3,1),cex=2)
  
  for(k in 1:length(cv_Folds))
  {
    # set.seed(201711800)
    k=3
    test_x = x[cv_Folds[[k]],]
    train_x = x[-cv_Folds[[k]],]
    
    test_y = y[cv_Folds[[k]],]
    train_y = y[-cv_Folds[[k]],]
    
    cvfit_cox = cv.glmnet(x, y, family = "cox",alpha=0.1,nlambda = 150)
    plot(cvfit_cox)
    
    ## train the model
    fit_cox = glmnet(train_x, train_y, family="cox",alpha=0.1,lambda = cvfit_cox$lambda.min)
    
    ## predict
    
    y_predict = predict(fit_cox,test_x,type = "response")
    
    
    data_a = data.frame(test_y,y_predict)
    colnames(data_a)=c('SurvObj','y_predict')
    
    # log rank & km curve
    data_a$signature = data_a$y_predict
    md_va = data_a$signature>1.5
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
    km_a = npsurv(SurvObj ~ signature, data = data_a)
    
    log_rank_predict_a=survdiff(SurvObj ~ signature, data = data_a)
    log_rank_predict_a_pvalue = pchisq(log_rank_predict_a$chisq, length(log_rank_predict_a$n)-1, lower.tail = FALSE)
    
    par(cex=1)
    survplot(km_a,xlab = 'months', n.risk   = TRUE,y.n.risk=0,cex.n.risk = 1,type="kaplan-meier",lty=1,lwd=3,
             label.curves = list(method = "arrow", cex = 1.2),col=c("orange","blue"),abbrev.label=FALSE)
    text(60,0.65,sprintf("p-value of log-rank test: %.10f", log_rank_predict_a_pvalue),col='brown2')
    title_K = paste("test:",k)
    title(title_K)
    
  }
  dev.off()
}
# dev.off()



