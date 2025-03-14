# 模型
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/Mime1包/mydata/GSE58812_data.Rdata")
Train$OS.time=as.numeric(Train$OS.time)
Train=Train %>% filter(OS.time>200)
#
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/自己代码跑/traindataset.Rdata")
TrainDataset$OS.time=as.numeric(TrainDataset$OS.time)
TrainDataset=TrainDataset %>% filter(OS.time>200)


sigGenes=c("OS.time","OS")
hub_gene=c(sel_genes,sigGenes)
Train=Train[,colnames(Train)%in%hub_gene]
data=Train
data$OS.time=as.numeric(data$OS.time)
data$OS=as.numeric(data$OS)
result=data.frame()
for(i in colnames(data[,3:ncol(data)])){
  cox<- coxph(Surv(OS.time, OS) ~ get(i), data = data)
  coxSummary = summary(cox)
  result=rbind(result,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}
result[,2:5] <-apply(result[,2:5],2,as.numeric)

# 通过P值以及HR对有预后意义的基因进行筛选

table(result$pvalue<0.05)
# FALSE  TRUE 
# 235    24
# outTab=result[result$pvalue<0.01,]
outTab=result[result$pvalue<0.05,]
genes=outTab$id
genes <- c("OS.time", "OS", genes)
#
# save(genes,file = "../单因素回归/手动单因素结果.Rdata")

# 森林图绘制：#读取输入文件
head(outTab)

data=outTab
# 将第一列变成行名：
rownames(data) <- data$id
data <- data[,-1]
head(data)

data$HR=as.numeric(data$HR)
data$HR.95L=as.numeric(data$HR.95L)
data$HR.95H=as.numeric(data$HR.95H)
data$pvalue=as.numeric(data$pvalue)
# rt <- data
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

gene <- rownames(rt)

n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse((as.numeric(hr)>1)&(pVal<0.05),'red3',"green3")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)

#
# 提取相同基因
comgene=Reduce(intersect,list(colnames(Train),colnames(TrainDataset),colnames(testDatasetGeo135565)))
Train=Train[,comgene]
TrainDataset=TrainDataset[,comgene]
testDatasetGeo135565=testDatasetGeo135565[,comgene]

# 只保留单因素回归后的基因
Train=Train[,genes]
TrainDataset=TrainDataset[,genes]
testDatasetGeo135565=testDatasetGeo135565[,genes]
#
str(Train)
Train$OS.time=as.numeric(Train$OS.time)
Train$OS=as.numeric(Train$OS)
#

TrainDataset$OS.time=as.numeric(TrainDataset$OS.time)
TrainDataset$OS=as.numeric(TrainDataset$OS)
testDatasetGeo135565$OS.time=as.numeric(testDatasetGeo135565$OS.time)
testDatasetGeo135565$OS=as.numeric(testDatasetGeo135565$OS)
#
trainlist=list(Train=Train,Test=TrainDataset
               # ,Test3=testDatasetGeo135565
)

#
result <- data.frame()
rf_nodesize <- 5
seed <- 123
#### 1-1.RSF #################
#################################################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = trainlist$Train,
             ntree = 1000,nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = trainlist$Train,
             ntree = best,nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)

rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)}) # 要改， list中一共有几个数据集就写几
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]

cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result,cc)

#################################################################
#### 1-2.rsf+Enet ################
#################################################################
#基因重要性
vi <- data.frame(imp=vimp.rfsrc(fit)$importance)
vi$imp <- (vi$imp-min(vi$imp))/(max(vi$imp)-min(vi$imp))
vi$ID <- rownames(vi)

#基因重要性可视化
pdf("rsf_highgene(rsf).pdf")
ggplot(vi,aes(imp,reorder(ID,imp)))+
  geom_bar(stat = 'identity',fill='#FF9933',color='black',width=0.7)+
  geom_vline(xintercept = 0.01,color='grey50',linetype=2)+
  labs(x='Relative importance by Random Forest',y=NULL)+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.x = element_text(size = 11,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.title = element_text(size=13,color='black'),
        legend.text = element_text(size=12,color='black'),
        legend.title = element_text(size=13,color='black'))+
  scale_y_discrete(expand = c(0.03,0.03))+
  scale_x_continuous(expand = c(0.01,0.01))
dev.off()
#提取重要性大于0.01的基因
rid <- rownames(vi)[vi$imp>0.01]
train2 <- Train[,c('OS.time','OS',rid)] # 这边记得要改 训练集叫啥，就写啥
trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})

x1 <- as.matrix(train2[,rid])
x2 <- as.matrix(Surv(train2$OS.time,train2$OS))

#利用循环探索Enet最佳模型
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))}) # 要改
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}
#利用交叉验证探索Enet最佳模型
set.seed(seed)
modelexp=as.matrix(trainlist2$Train[,c(3:ncol(trainlist2$Train))])
modelstat=Surv(trainlist2$Train$OS.time,trainlist2$Train$OS)
Enetmodel <- glmnet(modelexp,modelstat,family = 'cox',nfolds=10)
Enetmodel_cv<-cv.glmnet(modelexp,modelstat,family = 'cox',nfolds=10)
#建立最优模型
model_fit<-glmnet(modelexp,modelstat,family = 'cox',nfolds=10,keep=T,lambda = Enetmodel_cv$lambda.min)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model_fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=model_fit$lambda.min)))}) # 要改
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + Enet','[lambda=',round(Enetmodel_cv$lambda.min,3),']')
result <- rbind(result,cc)

####################################################################
####### 1-3.rsf+stepcox ####### #####
####################################################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,train2),direction = direction)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))}) # 要改
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + StepCox','[',direction,']')
  result <- rbind(result,cc)
}

#################################################################
#### 1-4.rsf+CoxBoost #### #######
#################################################################
set.seed(seed)
#计算最佳penalty
pen <- optimCoxBoostPenalty(train2[,'OS.time'],train2[,'OS'],as.matrix(train2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
#计算最佳stepno
cv.res <- cv.CoxBoost(train2[,'OS.time'],train2[,'OS'],as.matrix(train2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
#构建最佳模型
fit <- CoxBoost(train2[,'OS.time'],train2[,'OS'],as.matrix(train2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))}) # 要改
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result,cc)

#################################################################
#### 1-5.rsf+plsRcox #############
#################################################################
set.seed(seed)
#建立专用矩阵
model_exp=data.frame(train2[,-c(1:2)])
model_time=train2$OS.time
model_stat=train2$OS
#建立模型
model<-plsRcox(model_exp,time = model_time,event = model_stat,nt=10)
#进行交叉验证
cv.model<-cv.plsRcox(list(x=model_exp,time=model_time,status=model_stat),nt=5,verbose = T)
#构建最优模型
model<-plsRcox(model_exp,
               time = model_time,
               event = model_stat,
               nt=cv.model$lambda.min5,
               alpha.pvals.expli = 0.05,
               sparse = T,
               pvals.expli = T)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model,type="lp",newdata=x[,-c(1,2)])))}) # 要改
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + plsRcox')
result <- rbind(result,cc)

#################################################################
#### 1-6.rsf+superpc #############
#################################################################
data <- list(x=t(train2[,-c(1,2)]),y=train2$OS.time,censoring.status=train2$OS,featurenames=colnames(train2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) 
cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

rs <- lapply(trainlist2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + SuperPC')
result <- rbind(result,cc)

#################################################################
#### 1-7.rsf+gbm #################
#################################################################
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,
           data = train2,
           distribution = 'coxph',
           n.minobsinnode = 10,
           n.cores = 1,
           n.trees = 1000,
           shrinkage = 0.005,
           interaction.depth = 2,
           cv.folds = 5)
#构建最优模型
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,
           data = train2,
           distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,
           n.cores = 1)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))}) # 要改
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + GBM')
result <- rbind(result,cc)

#################################################################
#### 1-8.rsf+survivalsvm #########
#################################################################
set.seed(seed)
fit = survivalsvm(Surv(OS.time,OS)~., data= train2, gamma.mu = 2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))}) # 要改
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + survival-SVM')
result <- rbind(result,cc)

#################################################################
#### 1-9.rsf+Ridge ###############
#################################################################
set.seed(seed)
modelexp=as.matrix(train2[,c(3:ncol(train2))])
#利用循环探索Ridge最佳模型
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  model <- glmnet(modelexp,train2$OS,family = 'binomial',alpha = alpha,nfolds=10)
  model_cv<-cv.glmnet(modelexp,train2$OS,family = 'binomial',alpha =alpha,nfolds=10)
  fit<-glmnet(modelexp,train2$OS,family = 'binomial',alpha = alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="response",newx=as.matrix(x[,-c(1,2)]))))}) # 要改
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + Ridge','[α=',alpha,']')
  result <- rbind(result,cc)
}

# #################################################################
# #### 1-10.rsf+obliqueRSF #########
# #################################################################
# set.seed(seed)
# model<-orsf(data = train2,n_tree = 100,formula = Surv(OS.time,OS)~.)
# rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, new_data=x,pred_type = "risk")[,1]))}) # 要改
# rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
# rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
# rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
# rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
# cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
#   rownames_to_column('ID')
# cc$Model <- paste0('RSF + obliqueRSF')
# result <- rbind(result,cc)

#################################################################
#### 1-11.rsf+xgboost ############
#################################################################
set.seed(seed)
#建立专用矩阵
model_mat<-xgb.DMatrix(data = as.matrix(train2[,-c(1:2)]),label=train2$OS.time)
#构建参数
object<-list(bojective="surivival:cox",
             booster="gbtree",
             eval_metric="cox-nloglik",
             eta=0.01,
             max_depth=3,
             subsample=1,
             colsample_bytree=1,
             gamma=0.5)
#构建模型
model<-xgb.train(params=object,data = model_mat,nrounds = 100,watchlist = list(val2=model_mat),early_stopping_rounds = 10)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=as.matrix(x[,-c(1:2)]))))}) # 要改
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + xgboost')
result <- rbind(result,cc)

#################################################################
#### 1-12.rsf+CForest#############
#################################################################
set.seed(seed)
model<-party::cforest(Surv(OS.time,OS)~.,data=train2,controls = cforest_unbiased(ntree=50))
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CForest')
result <- rbind(result,cc)

#################################################################
#### 1-13.rsf+CTree###############
#################################################################
set.seed(seed)
#建立模型
model<-ctree(Surv(OS.time,OS)~.,data=train2)
#计算变量重要性
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CTree')
result <- rbind(result,cc)

#################################################################
#### 2-1.Enet ####################
#################################################################
modelexp=as.matrix(Train[,c(3:ncol(Train))]) # 要改，train
modelstat=Surv(Train$OS.time,Train$OS)
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  model <- glmnet(modelexp,modelstat,family = 'cox',alpha = alpha,nfolds=10)
  model_cv<-cv.glmnet(modelexp,modelstat,family = 'cox',alpha = alpha,nfolds=10)
  fit <-glmnet(modelexp,modelstat,family = 'cox',alpha =alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
  rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}

#################################################################
#### 2-2.Lasso+RSF################
#################################################################
set.seed(seed)
fit = cv.glmnet(modelexp, modelstat,family = "cox")
coef.min = coef(fit, s = "lambda.min") 
index=which(coef.min!=0)
actCoef=coef.min[index]
rid=row.names(coef.min)[index] # lasso筛选了
# 提取最佳模型的系数
best_model_coef <- coef(fit, s = fit$lambda.min)
# 10 x 1 sparse Matrix of class "dgCMatrix"
# 1
# BGN     0.05777262
# GADD45B .         
# JUN     0.41937235
# MMP14   .         
# NOTCH3  0.24966085
# PMEPA1  0.27412902
# PPDPF   .         
# SDC1    0.05284677
# TBX2    0.25064196
# UQCRQ   0.49472101
# rid <- coef.min@Dimnames[[1]] # 这取得是全部基因，没有筛选

train2 <- Train[,c('OS.time','OS',rid)] # 要改
trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})

set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = train2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = train2,
             ntree = best,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'Lasso + RSF'
result <- rbind(result,cc)

##################################################################
#### 2-3.Lasso+StepCox ############
##################################################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,train2),direction = direction)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + StepCox','[',direction,']')
  result <- rbind(result,cc)
}

####################################################################
#### 2-4.Lasso+CoxBoost #############
####################################################################
set.seed(seed)
#计算最佳penalty
modelpen<-optimCoxBoostPenalty(time = train2$OS.time,
                               status = train2$OS,
                               as.matrix(train2[,-c(1:2)]),
                               trace = T,
                               parallel = T)
#计算最佳stepno
cvmodel<-cv.CoxBoost(time = train2$OS.time,
                     status = train2$OS,
                     as.matrix(train2[,-c(1:2)]),
                     maxstepno = 100,
                     K = 3,
                     type = "verweij",
                     penalty=modelpen$penalty)
#构建CoxBoost模型
fit<-CoxBoost(time = train2$OS.time,
              status = train2$OS,
              as.matrix(train2[,-c(1:2)]),
              stepno = cvmodel$optimal.step,
              penalty = modelpen$penalty)

rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result,cc)

##################################################################
#### 2-5.Lasso+plsRcox ############
##################################################################
set.seed(seed)
#建立专用矩阵
model_exp=data.frame(train2[,-c(1:2)])
model_time=train2$OS.time
model_stat=train2$OS
#建立模型
model<-plsRcox(model_exp,time = model_time,event = model_stat,nt=10)
#进行交叉验证
cv.model<-cv.plsRcox(list(x=model_exp,time=model_time,status=model_stat),nt=5,verbose = F)
#构建最优模型
cv.plsRcox.res=cv.plsRcox(list(x=model_exp,time=model_time,status=model_stat),nt=5,verbose = F)
fit <- plsRcox(model_exp,
               time = model_time,
               event = model_stat,
               nt=cv.model$lambda.min5,
               alpha.pvals.expli = 0.05,
               sparse = T,
               pvals.expli = T)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + plsRcox')
result <- rbind(result,cc)

##################################################################
#### 2-6.Lasso+superpc ############
##################################################################
data <- list(x=t(train2[,-c(1,2)]),y=train2$OS.time,censoring.status=train2$OS,featurenames=colnames(train2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) 
cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(trainlist2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + SuperPC')
result <- rbind(result,cc)

#################################################################
#### 2-7.Lasso+gbm ###############
#################################################################
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = train2,distribution = 'coxph',
           n.trees = 1000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = train2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + GBM')
result <- rbind(result,cc)

#####################################################################
#### 2-8.Lasso+survivalsvm ###########
#####################################################################
set.seed(seed)
fit = survivalsvm(Surv(OS.time,OS)~., data= train2, gamma.mu = 2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + survival-SVM')
result <- rbind(result,cc)

#################################################################
#### 2-9.Lasso+Ridge #############
#################################################################
set.seed(seed)
modelexp=as.matrix(train2[,c(3:ncol(train2))])
#利用循环探索Ridge最佳模型
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  model <- glmnet(modelexp,train2$OS,family = 'binomial',alpha = alpha,nfolds=10)
  model_cv<-cv.glmnet(modelexp,train2$OS,family = 'binomial',alpha =alpha,nfolds=10)
  fit<-glmnet(modelexp,train2$OS,family = 'binomial',alpha = alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="response",newx=as.matrix(x[,-c(1,2)]))))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + Ridge','[α=',alpha,']')
  result <- rbind(result,cc)
}

# ###################################################################
# #### 2-10.Lasso+obliqueRSF #########
# ###################################################################
# set.seed(seed)
# model<-orsf(data = train2,n_tree = 100,formula = Surv(OS.time,OS)~.)
# rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, new_data=x,pred_type = "risk")[,1]))})
# rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
# rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
# rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
# rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
# cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
#   rownames_to_column('ID')
# cc$Model <- paste0('Lasso + obliqueRSF')
# result <- rbind(result,cc)

#################################################################
#### 2-11.Lasso+xgboost ##########
#################################################################
set.seed(seed)
#建立专用矩阵
model_mat<-xgb.DMatrix(data = as.matrix(train2[,-c(1:2)]),label=train2$OS.time)
#构建参数
object<-list(bojective="surivival:cox",
             booster="gbtree",
             eval_metric="cox-nloglik",
             eta=0.01,
             max_depth=3,
             subsample=1,
             colsample_bytree=1,
             gamma=0.5)
#构建模型
model<-xgb.train(params=object,data = model_mat,nrounds = 100,watchlist = list(val2=model_mat),early_stopping_rounds = 10)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=as.matrix(x[,-c(1:2)]))))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + xgboost')
result <- rbind(result,cc)

#################################################################
#### 2-12.Lasso+CForest###########
#################################################################
set.seed(seed)
model<-party::cforest(Surv(OS.time,OS)~.,data=train2,controls = cforest_unbiased(ntree=50))
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CForest')
result <- rbind(result,cc)

#################################################################
#### 2-13.Lasso+CTree#############
#################################################################
set.seed(seed)
model<-ctree(Surv(OS.time,OS)~.,data=train2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CTree')
result <- rbind(result,cc)

#################################################################
#### 3-1.StepCox #################
#################################################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction) # 要改
  rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}

#################################################################
#### 3-2.StepCox+RSF #############
#################################################################
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction) # 要改
  rid <- names(coef(fit))
  
  train2 <- Train[,c('OS.time','OS',rid)] # 要改
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~.,data = train2,
               ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  best <- which.min(fit$err.rate)
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~.,data = train2,
               ntree = best,nodesize = rf_nodesize,##该值建议多调整  
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + RSF')
  result <- rbind(result,cc)
}

#################################################################
#### 3-3.StepCox+Enet ###########
#################################################################
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  x1 <- as.matrix(train2[,rid])
  x2 <- as.matrix(Surv(train2$OS.time,train2$OS))
  
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
    rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
    rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
    rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
    rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
    rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox','[',direction,']',' + Enet','[α=',alpha,']')
    result <- rbind(result,cc)
  }
}
##################################################################
#### 3-4.StepCox+CoxBoost #########
##################################################################
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  
  set.seed(seed)
  pen <- optimCoxBoostPenalty(train2[,'OS.time'],train2[,'OS'],as.matrix(train2[,-c(1,2)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(train2[,'OS.time'],train2[,'OS'],as.matrix(train2[,-c(1,2)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  fit <- CoxBoost(train2[,'OS.time'],train2[,'OS'],as.matrix(train2[,-c(1,2)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + CoxBoost')
  result <- rbind(result,cc)
}

#################################################################
#### 3-5.StepCox+plsRcox #########
#################################################################
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  
  set.seed(seed)
  cv.plsRcox.res=cv.plsRcox(list(x=train2[,rid],time=train2$OS.time,status=train2$OS),nt=10,nfold = 10,verbose = F)
  fit <- plsRcox(train2[,rid],time=train2$OS.time,event=train2$OS,nt=as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + plsRcox')
  result <- rbind(result,cc)
}

#################################################################
#### 3-6.StepCox+superpc #########
#################################################################
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  
  data <- list(x=t(train2[,-c(1,2)]),y=train2$OS.time,censoring.status=train2$OS,featurenames=colnames(train2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) 
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                       n.fold = 5,
                       n.components=3,
                       min.features=1,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  rs <- lapply(trainlist2,function(w){
    test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2],RS=rr)
    return(rr2)
  })
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + SuperPC')
  result <- rbind(result,cc)
}

#############################################################
#### 3-7.StepCox+gbm #########
#############################################################
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = train2,distribution = 'coxph',
             n.trees = 1000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 5,n.cores = 1)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = train2,distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 5,n.cores = 1)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + GBM')
  result <- rbind(result,cc)
}

#######################################################################
#### 3-8.StepCox+survival-SVM ##########
#######################################################################
for (direction in c("both", "backward")) {
  #direction='both'
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  
  fit = survivalsvm(Surv(OS.time,OS)~., data= train2, gamma.mu = 1)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + survival-SVM')
  result <- rbind(result,cc)
}

################################################################
#### 3-9.StepCox+Ridge ##########
################################################################
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  set.seed(seed)
  modelexp=as.matrix(train2[,c(3:ncol(train2))])
  #利用循环探索Ridge最佳模型
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    model <- glmnet(modelexp,train2$OS,family = 'binomial',alpha = alpha,nfolds=10)
    model_cv<-cv.glmnet(modelexp,train2$OS,family = 'binomial',alpha =alpha,nfolds=10)
    fit<-glmnet(modelexp,train2$OS,family = 'binomial',alpha = alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
    rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="response",newx=as.matrix(x[,-c(1,2)]))))})
    rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
    rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
    rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
    rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox','[',direction,']','+ Ridge','[α=',alpha,']')
    result <- rbind(result,cc)
  }
}

# #######################################################################
# #### 3-10.StepCox+obliqueRSF ###########
# #######################################################################
# for (direction in c("both", "backward")) {
#   fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
#   rid <- names(coef(fit))
#   train2 <- Train[,c('OS.time','OS',rid)]
#   trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
#   model<-orsf(data = train2,n_tree = 100,formula = Surv(OS.time,OS)~.)
#   rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, new_data=x,pred_type = "risk")[,1]))})
#   rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
#   rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
#   rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
#   rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
#   cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
#     rownames_to_column('ID')
#   cc$Model <- paste0('StepCox','[',direction,']', '+ obliqueRSF')
#   result <- rbind(result,cc)
# }

###################################################################
#### 3-11.StepCox+xgboost ##########
###################################################################
for (direction in c("both", "backward")) {
  #direction='both'
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  
  #建立专用矩阵
  model_mat<-xgb.DMatrix(data = as.matrix(train2[,-c(1:2)]),label=train2$OS.time)
  #构建参数
  object<-list(bojective="surivival:cox",
               booster="gbtree",
               eval_metric="cox-nloglik",
               eta=0.01,
               max_depth=3,
               subsample=1,
               colsample_bytree=1,
               gamma=0.5)
  #构建模型
  model<-xgb.train(params=object,data = model_mat,nrounds = 100,watchlist = list(val2=model_mat),early_stopping_rounds = 10)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=as.matrix(x[,-c(1:2)]))))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']', '+ xgboost')
  result <- rbind(result,cc)
}

###################################################################
#### 3-12.StepCox+CForest ##########
###################################################################
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  model<-party::cforest(Surv(OS.time,OS)~.,data=train2,controls = cforest_unbiased(ntree=50))
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + CForest')
  result <- rbind(result,cc)
}

#################################################################
#### 3-13.StepCox+CTree ##########
#################################################################
for (direction in c("both", "backward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,Train),direction = direction)
  rid <- names(coef(fit))
  train2 <- Train[,c('OS.time','OS',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})
  model<-ctree(Surv(OS.time,OS)~.,data=train2)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=x,type = "response")))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + CTree')
  result <- rbind(result,cc)
}

############################################################
#### 4-1.CoxBoost ###########
############################################################
set.seed(seed)
pen <- optimCoxBoostPenalty(Train[,'OS.time'],Train[,'OS'],as.matrix(Train[,-c(1,2)]), # 要改
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(Train[,'OS.time'],Train[,'OS'],as.matrix(Train[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(Train[,'OS.time'],Train[,'OS'],as.matrix(Train[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)

###############################################################
#### 4-2.CoxBoost+Enet #########
###############################################################

rid <- names(coef(fit)[which(coef(fit)!=0)])
train2 <- Train[,c('OS.time','OS',rid)] # 要改
trainlist2 <- lapply(trainlist,function(x){x[,c('OS.time','OS',rid)]})

x1 <- as.matrix(train2[,rid])
x2 <- as.matrix(Surv(train2$OS.time,train2$OS))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}


##################################################################
#### 4-3.CoxBoost+stepcox #########
##################################################################

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,train2),direction = direction)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + StepCox','[',direction,']')
  result <- rbind(result,cc)
}

##############################################################
#### 4-4.CoxBoost+RSF #########
##############################################################

set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = train2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = train2,
             ntree = best,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'CoxBoost + RSF'
result <- rbind(result,cc)

#############################################################
#### 4-5.rsf+plsRcox #########
#############################################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train2[,rid],time=train2$OS.time,status=train2$OS),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(train2[,rid],time=train2$OS.time,event=train2$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + plsRcox')
result <- rbind(result,cc)

##################################################################
#### 4-6.CoxBoost+superpc #########
##################################################################
data <- list(x=t(train2[,-c(1,2)]),y=train2$OS.time,censoring.status=train2$OS,featurenames=colnames(train2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20, 
                     n.fold = 5,
                     n.components=3,
                     min.features=1,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(trainlist2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + SuperPC')
result <- rbind(result,cc)

##############################################################
#### 4-7.CoxBoost+gbm #########
##############################################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = train2,distribution = 'coxph',
           n.trees = 1000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 3,n.cores = 1)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = train2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 3,n.cores = 1)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + GBM')
result <- rbind(result,cc)

######################################################################
#### 4-8.CoxBoost+survivalsvm #########
######################################################################

fit = survivalsvm(Surv(OS.time,OS)~., data= train2, gamma.mu = 2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + survival-SVM')
result <- rbind(result,cc)

######################################################
#### 5.plsRcox#########
######################################################

set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=Train[,-c(1,2)],time=Train$OS.time,status=Train$OS),nt=10,nfold = 10,verbose = F) # 要改
fit <- plsRcox(Train[,-c(1,2)],time=Train$OS.time,event=Train$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)

######################################################
#### 6.superpc#########
######################################################

data <- list(x=t(Train[,-c(1,2)]),y=Train$OS.time,censoring.status=Train$OS,featurenames=colnames(Train)[-c(1,2)]) # 要改
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=1, 
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(trainlist,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result,cc)

###################################################
#### 7.GBM #########
###################################################

set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = Train,distribution = 'coxph', # 要改
           n.trees = 1000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = Train,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)

###########################################################
#### 8.survivalsvm #########
###########################################################

fit = survivalsvm(Surv(OS.time,OS)~., data= Train, gamma.mu = 2) # 要改
rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('survival-SVM')


result <- rbind(result,cc)
#
save(result,file = "../101机器学习/不用封装的跑/最终版_用于画图的result_24_9_11.Rdata")







##
rf_nodesize <- 5
seed <- 123
#
set.seed(seed)
#
Train[1:4,1:5]
#            OS.time OS      BGN   GADD45B      JUN
# GSM1419944    1066  1 10.40476  8.875197 11.23499
# GSM1419946     422  1 11.61945 10.102057 10.67574
# GSM1419948    4517  0 10.01371  9.349203 10.28792
# GSM1419949    3070  0 11.24955 10.199535 11.24470
training=Train
#
set.seed(seed)
# 1，确定最优penalty 
# 使用optimCoxBoostPenalty函数筛选当前数据的最优penalty ，将得到的pen$penalty 定为最终模型的参数
pen <- optimCoxBoostPenalty(training[,'OS.time'],
                            training[,'OS'],
                            as.matrix(training[,-c(1,2)]),
                            trace=TRUE,
                            start.penalty=500,
                            parallel = T)

pen$penalty
#[1] 5400
# 2，确定最优stepno
# 使用cv.CoxBoost函数确定最优stepno，取cv.res$optimal.step的值
#number of folds to be used for cross-validation
cv.res <- cv.CoxBoost(training[,'OS.time'],
                      training[,'OS'],
                      as.matrix(training[,-c(1,2)]),
                      maxstepno=500,
                      K=10,
                      type="verweij",
                      penalty= pen$penalty
                      # ,multicore=1
)

cv.res$optimal.step
#[1] 86
# 构建模型
# 使用上面得到的参数构建CoxBoost模型
fit <- CoxBoost(training[,'OS.time'],
                training[,'OS'],
                as.matrix(training[,-c(1,2)]),
                stepno=cv.res$optimal.step,
                penalty=pen$penalty)
summary(fit)
# 228 boosting steps resulting in 5 non-zero coefficients  
# partial log-likelihood: -70.00765 
# 
# parameter estimates > 0:
#   BGN, JUN, NOTCH3, PMEPA1, UQCRQ 
# parameter estimates < 0:
plot(fit)
selected_vars_cb <- names(coef(fit))[coef(fit) != 0]
# 查看选择的变量
print(selected_vars_cb)

data=training
# 运行Enet[0]
data_coxboost <- data %>% dplyr::select(c("OS.time", "OS", selected_vars_cb))
#
x <- as.matrix(data_coxboost[,c(3:ncol(data_coxboost))])   # 预测变量矩阵（除去 OS.time 和 OS 列）
y <- Surv(data_coxboost$OS.time,data_coxboost$OS)  # 生存对象
# 3. 使用 Elastic Net 

alpha_value <- 0.5
# 使用 glmnet 进行模型拟合
enet_model <- glmnet(x, y, family = "cox", alpha = alpha_value,nfolds=10)
cv_enet <- cv.glmnet(x, y, family = "cox", alpha = alpha_value,nfolds=10)
cv_enet
best_lambda <- cv_enet$lambda.min
plot(cv_enet)
# coef_enet <- coef(enet_model, s = best_lambda)
# coef_enet=as.matrix(coef_enet)
# selected_vars_enet <- rownames(coef_enet)[coef_enet != 0]
# print(selected_vars_enet)



#此处使用lambda.min, 也可以尝试lambda.1se
coefficient <- coef(cv_enet, s = best_lambda) 

#系数不等于0的为纳入的变量（基因）
Active.index <- which(as.numeric(coefficient) != 0)
Active.coefficient <- as.numeric(coefficient)[Active.index]
sig_gene_mult_cox <- rownames(coefficient)[Active.index]
#查看具体哪些基因
sig_gene_mult_cox
training_cox <- data_coxboost %>% 
  dplyr::select(OS,OS.time,sig_gene_mult_cox) 
multiCox <- coxph(Surv(OS.time, OS) ~ ., data =  training_cox)
coefficients <- coef(multiCox)
selected_vars_stepcox <- names(coefficients)
cat("Selected Variables after StepCox:", paste(selected_vars_stepcox, collapse = ", "), "\n")
print(coefficients)

#predict函数计算风险评分
riskScore=predict(multiCox,type="risk",newdata=training_cox) 
riskScore<-as.data.frame(riskScore)
riskScore$sample <- rownames(riskScore)
head(riskScore,2) 
#
riskScore_cli <- training_cox %>% 
  rownames_to_column("sample") %>% 
  inner_join(riskScore)
# 计算ROC
data=riskScore_cli
data$risk_score=data$riskScore
time_points <- c( 3, 5,7) * 365

roc_3 <- survivalROC(Stime = data$OS.time, status = data$OS, marker = data$risk_score, predict.time = time_points[1], method = "KM")
roc_5 <- survivalROC(Stime = data$OS.time, status = data$OS, marker = data$risk_score, predict.time = time_points[2], method = "KM")
roc_7 <- survivalROC(Stime = data$OS.time, status = data$OS, marker = data$risk_score, predict.time = time_points[3], method = "KM")

riskScore_cli$riskScore2 <- ifelse(riskScore_cli$riskScore > median(riskScore_cli$riskScore),
                                   "High","Low")
#KM分析
fit2 <- survfit(Surv(OS.time, as.numeric(OS)) ~ riskScore2, data=riskScore_cli)

p2 <- ggsurvplot(fit2, data = riskScore_cli,
                 pval = T,
                 risk.table = T,
                 surv.median.line = "hv", #添加中位生存曲线
                 palette=c("red", "blue"),  #更改线的颜色
                 legend.labs=c("High risk","Low risk"), #标签
                 legend.title="RiskScore",
                 title="Overall survival", #标题
                 ylab="Cumulative survival (percentage)",xlab = " Time (Days)", #更改横纵坐标
                 censor.shape = 124,censor.size = 2,conf.int = FALSE, #删失点的形状和大小 break.x.by = 720#横坐标间隔 
) 
p2
#
Train_riskScore_cli=riskScore_cli
# 验证集呢
traindataset_cox <- TrainDataset %>% 
  dplyr::select(OS,OS.time,sig_gene_mult_cox)
riskScore=predict(multiCox,type="risk",newdata=traindataset_cox) 
riskScore<-as.data.frame(riskScore)
riskScore$sample <- rownames(riskScore)
head(riskScore,2) 
#
riskScore_cli <- traindataset_cox %>% 
  rownames_to_column("sample") %>% 
  inner_join(riskScore)
data=riskScore_cli
data$risk_score=data$riskScore
time_points <- c(3, 5,7) * 365

roc_3 <- survivalROC(Stime = data$OS.time, status = data$OS, marker = data$risk_score, predict.time = time_points[1], method = "KM")
roc_5 <- survivalROC(Stime = data$OS.time, status = data$OS, marker = data$risk_score, predict.time = time_points[2], method = "KM")
roc_7 <- survivalROC(Stime = data$OS.time, status = data$OS, marker = data$risk_score, predict.time = time_points[3], method = "KM")

riskScore_cli$riskScore2 <- ifelse(riskScore_cli$riskScore > median(riskScore_cli$riskScore),
                                   "High","Low")
#KM分析
fit2 <- survfit(Surv(OS.time, as.numeric(OS)) ~ riskScore2, data=riskScore_cli)

p2 <- ggsurvplot(fit2, data = riskScore_cli,
                 pval = T,
                 risk.table = T,
                 surv.median.line = "hv", #添加中位生存曲线
                 palette=c("red", "blue"),  #更改线的颜色
                 legend.labs=c("High risk","Low risk"), #标签
                 legend.title="RiskScore",
                 title="Overall survival", #标题
                 ylab="Cumulative survival (percentage)",xlab = " Time (Days)", #更改横纵坐标
                 censor.shape = 124,censor.size = 2,conf.int = FALSE, #删失点的形状和大小 break.x.by = 720#横坐标间隔 
) 
p2

