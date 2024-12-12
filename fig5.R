####Fig5I####
{
  ##### 01-诺莫图 #####
  # 接下来是真实数据，然后使用普通画法
  # 需要修改, 全部改成因子类型
  riskScore_cli2$Stage <- ifelse(riskScore_cli2$Stage  %in%   c("I","II"),  "Low", "High")
  riskScore_cli2$Tumor <- ifelse(riskScore_cli2$Tumor  %in%   c("T1","T4"), "High", "Low")
  
  
  riskScore_cli2 <- Test_riskScore_cli_4
  # 设置因子型变量
  riskScore_cli2$Stage <- factor(riskScore_cli2$Stage, levels = c("High","Low" ))
  levels(riskScore_cli2$Stage)
  
  riskScore_cli2$Tumor <- factor(riskScore_cli2$Tumor, levels = c("High","Low" ))
  levels(riskScore_cli2$Tumor)
  
  
  riskScore_cli2$Metastasis <- factor(riskScore_cli2$Metastasis, levels = c("Yes","No"))
  levels(riskScore_cli2$Metastasis)
  
  riskScore_cli2$Weight <- factor(riskScore_cli2$Weight, levels = c("High","Low" ))
  levels(riskScore_cli2$Weight)
  
  
  riskScore_cli2$riskScore2 <- factor(riskScore_cli2$riskScore2, levels = c("High","Low" ))
  levels(riskScore_cli2$riskScore2)
  
  
  riskScore_cli2$Age <- factor(riskScore_cli2$Age, levels = c("High","Low" ))
  
  
  levels(riskScore_cli2$Age)
  
  
  
  
  
  colnames(riskScore_cli2)[c(5,6)] <- c( "OS", "OS.time")
  riskScore_cli2$OS.time <- round(riskScore_cli2$OS.time/365, 2)
  
  colnames(riskScore_cli2)
  
  # 然后发现没有临床信息，所以需要加载GSEA的信息
  dd <- datadist(riskScore_cli2)
  options(datadist="dd")
  
  
  f <- psm(Surv(OS.time,OS) ~ Stage + Tumor +  Weight + riskScore2 + Metastasis+age,
           data = riskScore_cli2,dist='lognormal')
  surv <- Survival(f) # 构建生存概率函数
  ## time是以”天“为单位,此处绘制1年，3,5年的生存概率
  nom <- nomogram(f, fun=list(function(x) surv(1, x),
                              function(x) surv(3, x),
                              function(x) surv(5, x) ),
                  funlabel=c("1-year OS", "3-year OS", 
                             "5-year OS"))
  plot(nom, xfrac=.2)
  
  
  
  # 接下来使用优美画法
  
  library(regplot)
  Cox_nomo2 <- cph(Surv(OS.time,OS) ~ Stage + Tumor +  Weight + riskScore2 + Metastasis+age ,
                   data = riskScore_cli2,dist='lognormal', x=T, y=T)
  
  regplot(Cox_nomo2, 
          # observation = riskScore_cli2[4,], #指定某一患者
          interval ="confidence", 
          title="Nomogram",
          plots=c("violin", "boxes"), 
          clickable = T,
          failtime = c(1,3,5)) #设置随访时间1年、3年和5年
  
}

#### Fig5GH####
{
  ##### 02-COX曲线 #####
  cox_need
  library(tidyverse)
  cox_need <- rownames_to_column(cox_need, var = "Variable")
  cox_need$hazard_ratio <- paste0(cox_need$HR,"(",cox_need$ower_95,"-",cox_need$upper_95,")")
  
  dat <- cox_need[,c(1,6,10,11,12,5)]
  colnames(dat)=c("Variable","HR","CI5","CI95","HR (95% CI)","Pvalue")
  dat = rbind(c("Names", NA,NA,NA,"Hazard Ratio(95% CI)", "p.value"),dat)
  dat <- dat[c(1,7,6,4,5,3,2),]
  
  
  
  
  dat$Variable <- c(NA,"Age","Weight","Stage","Tumor","Metastasis","RiskScore")
  dat <- dat[,c(1:4,6,5)]
  
  #画图
  forestplot(dat[,c(1,5,6)], #显示表格的第1，5，6列内容
             mean=dat[,2],   #第2列为HR，变成森林图的方块
             lower=dat[,3], upper=dat[,4], #第3列为5%CI，第4列为95%CI，将化作线段，穿过方块
             zero=1,            #零线或参考线位置为HR=1
             boxsize=0.2,       #设置方块大小
             graph.pos=4,#将森林图插在第3列
             xticks=c(-5,0,5,10,15) ,# 设置横轴数字
             txt_gp=fpTxtGp (
               label=gpar(cex=0.8) ,ticks=gpar(cex=0.6)
             ),#调整字体
             # hrzl_lines=list("1" = gpar(lty=1,lwd=1.5),
             #                 "2" = gpar(lty=1,lwd=1.5),
             #                 "5"= gpar(lty=1,lwd=1.5)), # 在1,2,7行添加横线
             col=fpColors ( box = 'red ' , #方块颜色 
                            lines = ' black ' ,#置信区间横线颜色
                            zero = "grey" ),#参考线颜色
             lwd.zero=1,#参考线宽度
             lwd.ci=1.5, # 置信区间横线宽度
             lty.ci=7 ,# 置信区间横线类型
             ci.vertices.height=0.1,   #置信区间横线两端竖线高度  
             title = "OS (Univariate Cox)"
  )
  
  
  
  # 单因素出完出多因素的
  # coxph(Surv(futime, fustat) ~ Metastasis + RiskScore + Stage, data = Test_riskScore_cli_4)
  colnames(Test_riskScore_cli_4)
  res <- coxph(Surv(futime, fustat) ~ Metastasis + riskScore2 + Stage, data = Test_riskScore_cli_4)
  summary(res)
  
  mul_cox <- summary(res)
  mul_HR<- round(mul_cox$coefficients[,2],2) 
  mul_Pvalue<- round(mul_cox$coefficients[,5],4) 
  mul_CI5<-round(mul_cox$conf.int[,3],2)
  mul_CI95<-round(mul_cox$conf.int[,4],2)
  mul_CI<-paste0(mul_HR,' (',mul_CI5,'-',mul_CI95,')')
  Variable<-row.names(data.frame(mul_cox$coefficients))
  mulcox_res<- data.frame(Variable,mul_HR,mul_CI5,mul_CI95,mul_CI,mul_Pvalue)
  colnames(mulcox_res)=c("Variable","HR","CI5","CI95","HR (95% CI)","Pvalue")
  View(mulcox_res)
  
  library(forestplot)
  #添加表头
  dat=rbind(c("Variable", NA,NA,NA,"HR (95% CI)", "Pvalue"),mulcox_res)
  dat <- dat[,c(1:4,6,5)]
  
  dat[1,1] <- NA
  dat[2,1] <- "Metastasis"
  dat[3,1] <- "RiskScore"
  dat[4,1] <- "Stage"
  #画图
  forestplot(dat[,c(1,5,6)], #显示表格的第1，5，6列内容
             mean=dat[,2],   #第2列为HR，变成森林图的方块
             lower=dat[,3], upper=dat[,4], #第3列为5%CI，第4列为95%CI，将化作线段，穿过方块
             zero=1,            #零线或参考线位置为HR=1
             boxsize=0.2,       #设置方块大小
             graph.pos=4,#将森林图插在第3列
             xticks=c(-5,0,5,10,15) ,# 设置横轴数字
             txt_gp=fpTxtGp (
               label=gpar(cex=0.8) ,ticks=gpar(cex=0.6)
             ),#调整字体
             # hrzl_lines=list("1" = gpar(lty=1,lwd=1.5),
             #                 "2" = gpar(lty=1,lwd=1.5),
             #                 "5"= gpar(lty=1,lwd=1.5)), # 在1,2,7行添加横线
             col=fpColors ( box = 'red ' , #方块颜色 
                            lines = ' black ' ,#置信区间横线颜色
                            zero = "grey" ),#参考线颜色
             lwd.zero=1,#参考线宽度
             lwd.ci=1.5, # 置信区间横线宽度
             lty.ci=7 ,# 置信区间横线类型
             ci.vertices.height=0.1,   #置信区间横线两端竖线高度  
             title = "OS (Multivariate Cox)"
  )
  
}

####Fig5J####
{
  #####03-timeROC#####
  install.packages('timeROC')
  library(timeROC)
  
  time_points <- c(1, 3, 5, 7) * 365.25  # 将年份转换为天数
  
  data <- Test_riskScore_cli_4
  data$OS.time <- data$futime
  data$OS <- data$fustat
  data$risk_score <- data$riskScore
  # 计算 1 年、3 年、5 年的 ROC 曲线
  roc_1 <- timeROC(T = data$OS.time, delta = data$OS, marker = data$risk_score, cause = 1, times = time_points[1])
  roc_3 <- timeROC(T = data$OS.time, delta = data$OS, marker = data$risk_score, cause = 1, times = time_points[2])
  roc_5 <- timeROC(T = data$OS.time, delta = data$OS, marker = data$risk_score, cause = 1, times = time_points[3])
  roc_7 <- timeROC(T = data$OS.time, delta = data$OS, marker = data$risk_score, cause = 1, times = time_points[4])
  roc_1$AUC;roc_3$AUC;roc_5$AUC;roc_7$AUC
  
  data$fustat <- as.factor(data$fustat)
  data$futime<- round(data$futime/365, 2)
  ## 构建timeROC
  ROC <- timeROC(T=data$futime, #生存时间
                 delta=data$fustat,   #生存状态
                 marker=data$riskScore, #计算timeROC的变量
                 cause=1,                #阳性结局指标数值(1表示死亡)
                 weighting="marginal",   #计算方法，默认为marginal
                 times=c(1, 3, 5, 7),       #时间点，选取1年，3年和5年的生存率
                 iid=TRUE)
  ROC
  plot(ROC,
       time=1, col="red", lty=1,lwd=2, title = "")   #time是时间点，col是线条颜色、lty为图例线条类型、lwd为图例线条宽度
  plot(ROC,
       time=3, col="blue", add=TRUE, lty=1,lwd=2)    #add指是否添加在上一张图中
  plot(ROC,
       time=5, col="orange", add=TRUE, lty=1,lwd=2)
  plot(ROC,
       time=7, col="yellow", add=TRUE, lty=1,lwd=2)
  ## 添加图例
  legend("bottomright",#图例画在右下角
         c(paste0("AUC at 1 year: ",round(ROC[["AUC_1"]][1],2)), #提取1年AUC构建图例标签
           paste0("AUC at 3 year: ",round(ROC[["AUC_1"]][2],2)), #提取3年AUC构建图例标签
           paste0("AUC at 5 year: ",round(ROC[["AUC_1"]][3],2)),#提取5年AUC构建图例标签
           paste0("AUC at 7 year: ",round(ROC[["AUC_1"]][4],2))),#提取5年AUC构建图例标签
         col=c("red",
               "blue",
               "orange",
               "yellow"), #设置1，3，5年AUC图例标签的图例颜色，注意与曲线保持对应
         lty=1,  
         lwd=2,  
         bty = "n" #o表示用框框把图例部分框起来，为默认。n表示不画框框
  )
  
  
  
  colnames(Test_riskScore_cli_4)
  f <- coxph(Surv(futime, fustat)~riskScore + Stage +  Metastasis, data = Test_riskScore_cli_4)
  
  Test_riskScore_cli_4$riskscoreall <- predict(f, type ="lp")
  
  
  
  data <- Test_riskScore_cli_4
  data$fustat <- as.factor(data$fustat)
  
  ## 构建timeROC
  ROC <- timeROC(T=data$futime, #生存时间
                 delta=data$fustat,   #生存状态
                 marker=data$riskscoreall, #计算timeROC的变量
                 cause=1,                #阳性结局指标数值(1表示死亡)
                 weighting="marginal",   #计算方法，默认为marginal
                 times=c(1, 3, 5, 7)*365,       #时间点，选取1年，3年和5年的生存率
                 iid=TRUE)
  ROC
  plot(ROC,
       time=1*365, col="red", lty=1,lwd=2, title = "")   #time是时间点，col是线条颜色、lty为图例线条类型、lwd为图例线条宽度
  plot(ROC,
       time=3*365, col="blue", add=TRUE, lty=1,lwd=2)    #add指是否添加在上一张图中
  plot(ROC,
       time=5*365, col="orange", add=TRUE, lty=1,lwd=2)
  plot(ROC,
       time=7*365, col="yellow", add=TRUE, lty=1,lwd=2)
  ## 添加图例
  legend("bottomright",#图例画在右下角
         c(paste0("AUC at 1 year: ",round(ROC[["AUC_1"]][1],2)), #提取1年AUC构建图例标签
           paste0("AUC at 3 year: ",round(ROC[["AUC_1"]][2],2)), #提取3年AUC构建图例标签
           paste0("AUC at 5 year: ",round(ROC[["AUC_1"]][3],2)),#提取5年AUC构建图例标签
           paste0("AUC at 7 year: ",round(ROC[["AUC_1"]][4],2))),#提取5年AUC构建图例标签
         col=c("red",
               "blue",
               "orange",
               "yellow"), #设置1，3，5年AUC图例标签的图例颜色，注意与曲线保持对应
         lty=1,  
         lwd=2,  
         bty = "n" #o表示用框框把图例部分框起来，为默认。n表示不画框框
  )
  

  library(survival)
  data$facriskscoreall <- ifelse(data$riskscoreall>=median(data$riskscoreall), "High", "Low")
  data$facriskscoreall <- factor(data$facriskscoreall, levels = c("Low", "High"))
  
  # cox <- coxph(Surv(futime,fustat)~riskScore2 + Metastasis + Stage, data = data)
  # summary(cox)
  
  
  cox <- coxph(Surv(futime, fustat)~riskScore2 + Stage +  Metastasis, data = Test_riskScore_cli_4)
  summary(cox)
  
  # 然后基于riskscore和临床指标预测临床风险评分：
  data$riskscore3 <- predict(cox, type = "lp")
  head(data)
  # riskscore2 OS.time OS riskscore3
  # 1   1.634682     967  1   1.634682
  # 2   1.634682     584  0   1.634682
  # 3   1.634682    2654  0   1.634682
  # 4   0.000000     754  1   0.000000
  # 5   0.000000    2048  0   0.000000
  # 6   1.634682    1027  0   1.634682
  
  
  
  
  
  # 然后进行ROC分析：
  
  
  
  # 
  # data$OS.time <- data$futime
  # data$OS <- data$fustat
  
  library(survivalROC)
  
  ROC1<- survivalROC(Stime=data$OS.time, 
                     status=data$OS, 
                     marker = data$riskscore3, 
                     predict.time =365*1, #1年生存率,也可改为3年生存率 365*3
                     method = "KM")
  
  ROC1
  
  
  
  
  
  
  ROC2<- survivalROC(Stime=data$OS.time, 
                     status=data$OS, 
                     marker = data$riskscore3, 
                     predict.time =365*3, #1年生存率,也可改为3年生存率 365*3
                     method = "KM")
  
  ROC2
  
  
  
  ###############画ROC3
  ROC3<- survivalROC(Stime=data$OS.time, 
                     status=data$OS, 
                     marker = data$riskscore3, 
                     predict.time =365*5, #1年生存率,也可改为3年生存率 365*3
                     method = "KM")
  
  ROC3
  
  
  
  
  
  
  #
  ROC4<- survivalROC(Stime=data$OS.time, 
                     status=data$OS, 
                     marker = data$riskscore3, 
                     predict.time =365*7, #1年生存率,也可改为3年生存率 365*3
                     method = "KM")
  
  ROC4
  
  # 画图所有代码在这里
  {
    
    plot(ROC1$FP, ROC1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red", 
         xlab="False positive rate", ylab="True positive rate",
         lwd = 2, cex.main=1, cex.lab=1, cex.axis=1.2, font=1.2)
    abline(0,1)
    aucText=paste0("1 years"," (AUC=",sprintf("%.3f",ROC1$AUC),")") #这个后面添加legend用
    
    
    lines(ROC2$FP, ROC2$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="orange",lwd = 2)
    aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",ROC2$AUC),")") #这个后面添加legend用
    
    lines(ROC3$FP, ROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="cyan",lwd = 2)
    aucText5=paste0("5 years"," (AUC=",sprintf("%.3f",ROC3$AUC),")") #这个后面添加legend用
    
    lines(ROC4$FP, ROC4$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="yellow",lwd = 2)
    aucText7=paste0("7 years"," (AUC=",sprintf("%.3f",ROC4$AUC),")") #这个后面添加legend用
    #添加legend
    legend("bottomright", c(aucText,aucText3,aucText5,aucText7),
           lwd=2,bty="n",col=c("red","orange","cyan","yellow"),cex=0.8)
  }
  
  
  
  
  ROC1$AUC;ROC2$AUC;ROC3$AUC;ROC4$AUC
  
  
  plot(ROC1,
       time=1*365, col="red", lty=1,lwd=2, title = "")   #time是时间点，col是线条颜色、lty为图例线条类型、lwd为图例线条宽度
  plot(ROC2,
       time=3*365, col="blue", add=TRUE, lty=1,lwd=2)    #add指是否添加在上一张图中
  plot(ROC3,
       time=5*365, col="orange", add=TRUE, lty=1,lwd=2)
  plot(ROC4,
       time=7*365, col="yellow", add=TRUE, lty=1,lwd=2)
  ## 添加图例
  legend("bottomright",#图例画在右下角
         c(paste0("AUC at 1 year: ",round(ROC[["AUC_1"]][1],2)), #提取1年AUC构建图例标签
           paste0("AUC at 3 year: ",round(ROC[["AUC_1"]][2],2)), #提取3年AUC构建图例标签
           paste0("AUC at 5 year: ",round(ROC[["AUC_1"]][3],2)),#提取5年AUC构建图例标签
           paste0("AUC at 7 year: ",round(ROC[["AUC_1"]][4],2))),#提取5年AUC构建图例标签
         col=c("red",
               "blue",
               "orange",
               "yellow"), #设置1，3，5年AUC图例标签的图例颜色，注意与曲线保持对应
         lty=1,  
         lwd=2,  
         bty = "n" #o表示用框框把图例部分框起来，为默认。n表示不画框框
  )
}

####Fig5K####
{
  #####04-C-index#####
  # C-index
  library(rms)
  library(pec)
  library(ggplot2)
  riskScore_cli2$OS.time <- riskScore_cli2$OS.time2
  
  
  ###
  models=list(  Stage=cph(Surv(OS.time,OS)~Stage,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                Tumor=cph(Surv(OS.time,OS)~Tumor,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                Weight=cph(Surv(OS.time,OS)~Weight,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                riskScore2=cph(Surv(OS.time,OS)~riskScore2,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                Metastasis=cph(Surv(OS.time,OS)~Metastasis,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                Age=cph(Surv(OS.time,OS)~Age,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                riskscoreall=cph(Surv(OS.time,OS)~riskscoreall,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T))
  
  
  #time-cindex计算
  times<-c(1,3,5,7)
  cindex<-cindex(models,
                 formula=Surv(OS.time,OS)~1,
                 eval.times=times,
                 data=riskScore_cli2)
  
  plot(cindex)
  
  
  cindex$AppCindex
  
  
  cindex_df<-data.frame(
    Time=times,
    do.call(cbind,cindex$AppCindex)
  )
  cindex_df
  
  library(tidyr)
  dat=pivot_longer(cindex_df,cols=2:7,
                   names_to="model",
                   values_to="cindex")
  head(dat)
  
  
  library(ggplot2)
  ggplot(dat,aes(x=Time,y=cindex))+
    geom_line(aes(color=model),linewidth=2)+
    scale_color_brewer(palette="Set1")+
    geom_hline(yintercept=0.5,linetype=4)+
    ylim(0.4,1)+
    labs(title="Time-dependentC-index",x="Time(years)",y="C-index")+
    theme_bw()
  
  
  
  
  

  library(survival)
  
  data <- riskScore_cli2
  cox <- coxph(Surv(OS.time,OS)~as.factor(riskScore2) + as.factor(Stage),data)
  summary(cox)
  
  
  # 然后基于riskscore和临床指标预测临床风险评分：
  data$riskscoreall <- predict(cox, type = "lp")
  head(data)
  
  data$riskscoreall <- ifelse(data$riskscoreall >= median(data$riskscoreall), "High", "Low")
  data$riskscoreall <- factor(data$riskscoreall, levels = c("Low","High"))
  
  
  riskScore_cli2 <- data
  models=list(  Tumor=cph(Surv(OS.time,OS)~Tumor,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                Weight=cph(Surv(OS.time,OS)~Weight,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                riskScore2=cph(Surv(OS.time,OS)~riskScore2,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                Metastasis=cph(Surv(OS.time,OS)~Metastasis,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                Age=cph(Surv(OS.time,OS)~Age,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T),
                riskscoreall=cph(Surv(OS.time,OS)~riskscoreall,data=riskScore_cli2,x=TRUE,y=TRUE,surv=T))
  
  #time-cindex计算
  times<-c(1,3,5,7)
  cindex<-cindex(models,
                 formula=Surv(OS.time,OS)~1,
                 eval.times=times,
                 data=riskScore_cli2)
  
  plot(cindex)
  
  
  cindex$AppCindex
  
  
  cindex_df<-data.frame(
    Time=times,
    do.call(cbind,cindex$AppCindex)
  )
  cindex_df
  
  library(tidyr)
  dat=pivot_longer(cindex_df,cols=2:7,
                   names_to="model",
                   values_to="cindex")
  head(dat)
  
  
  library(ggplot2)
  ggplot(dat,aes(x=Time,y=cindex))+
    geom_line(aes(color=model),linewidth=2)+
    scale_color_brewer(palette="Set1")+
    geom_hline(yintercept=0.5,linetype=4)+
    ylim(0.4,1)+
    labs(title="Time-dependentC-index",x="Time(years)",y="C-index")+
    theme_bw()
  
}

####Fig5L####
{
  #####05-决策曲线 ######
  library(rms)
  Stage <- lrm(OS~Stage,riskScore_cli2)
  Metastasis <- lrm(OS~Metastasis,riskScore_cli2)
  Nomogram <- lrm(OS~Stage+Metastasis+riskScore2,riskScore_cli2)
  
  d_train <- dca(Stage,Metastasis,Nomogram)
  ggplot(d_train)
}
