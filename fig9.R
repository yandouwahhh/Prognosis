# cellchat same as fig3
# KM
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/最相关模块基因.Rdata")
# 过滤样本
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/Mime1包/mydata/GSE58812_不分割log的清洗好的data.Rdata")
Train$OS.time=as.numeric(Train$OS.time)
Train=Train %>% filter(OS.time>200)
#
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/自己代码跑/traindataset.Rdata")
TrainDataset$OS.time=as.numeric(TrainDataset$OS.time)
TrainDataset=TrainDataset %>% filter(OS.time>200)
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/自己代码跑/testDatasetGeo135565.Rdata")
testDatasetGeo135565$OS.time=as.numeric(testDatasetGeo135565$OS.time)
testDatasetGeo135565=testDatasetGeo135565 %>% filter(OS.time>200)
#
Train$OS.time=as.numeric(Train$OS.time)
Train$OS=as.numeric(Train$OS)
#

TrainDataset$OS.time=as.numeric(TrainDataset$OS.time)
TrainDataset$OS=as.numeric(TrainDataset$OS)
testDatasetGeo135565$OS.time=as.numeric(testDatasetGeo135565$OS.time)
testDatasetGeo135565$OS=as.numeric(testDatasetGeo135565$OS)


library(survival)
library(survminer)
#
Train[1:4,1:4]
#                    ID OS.time OS     A1BG
# GSM1419942 GSM1419942    1520  1 8.591714
# GSM1419943 GSM1419943    1281  1 7.533359
# GSM1419944 GSM1419944    1066  1 7.366889
# GSM1419945 GSM1419945    1050  1 7.574151
exp_meta_for_survival=Train
# 按基因表达值，计算二分类
exp_meta_for_survival$SULF1_Group=ifelse(exp_meta_for_survival$WFDC1>median(exp_meta_for_survival$WFDC1),'high','low')
exp_meta_for_survival$SULF1_Group=ifelse(exp_meta_for_survival$UQCRQ>median(exp_meta_for_survival$UQCRQ),'high','low')
exp_meta_for_survival$SULF1_Group=ifelse(exp_meta_for_survival$ISCU>median(exp_meta_for_survival$ISCU),'high','low')
exp_meta_for_survival$SULF1_Group=ifelse(exp_meta_for_survival$FN1>median(exp_meta_for_survival$FN1),'high','low')
exp_meta_for_survival$SULF1_Group=ifelse(exp_meta_for_survival$TBX2>median(exp_meta_for_survival$TBX2),'high','low')
exp_meta_for_survival$SULF1_Group=ifelse(exp_meta_for_survival$NOTCH3>median(exp_meta_for_survival$NOTCH3),'high','low')
exp_meta_for_survival$SULF1_Group=ifelse(exp_meta_for_survival$C1S>median(exp_meta_for_survival$C1S),'high','low')
exp_meta_for_survival$SULF1_Group=ifelse(exp_meta_for_survival$JUN>median(exp_meta_for_survival$JUN),'high','low')
exp_meta_for_survival$SULF1_Group=ifelse(exp_meta_for_survival$PLS3>median(exp_meta_for_survival$PLS3),'high','low')

table(exp_meta_for_survival$SULF1_Group)

# sfit <- survfit(Surv(OS.time, OS)~SULF1_Group, data=exp_meta_for_survival)
# print(sfit)
# ggsurvplot(sfit, conf.int=F, pval=TRUE)
#
fit <- survfit(Surv(OS.time, OS) ~ SULF1_Group, data = exp_meta_for_survival)
ggsurvplot(fit, data = exp_meta_for_survival,  
           conf.int = TRUE,  pval = TRUE, 
           surv.median.line = "hv", 
           risk.table = TRUE, palette = "hue") 

library(ggsci)
palette = "npg"  # 使用 Nature 期刊风格配色

ggsurvplot(fit, data = exp_meta_for_survival,
           conf.int = TRUE,  pval = TRUE,
           surv.median.line = "hv",
           # risk.table = TRUE,
           # risk.table.height = 0.25,  # 调整风险表高度
           # risk.table.col = "strata", # 风险表按分层显示颜色
           palette = "npg",
           legend.labs=c("High","Low"), #标签
           legend.title="RiskScore",
           title="PLS3", #标题
           ylab="Overall survival",xlab = " Time (Days)", #更改横纵坐标
           censor.shape = 124,censor.size = 2, #删失点的形状和大小 break.x.by = 720#横坐标间隔 
           ggtheme = theme_minimal()+ # 使用干净的主题
             theme(plot.title = element_text(hjust = 0.5),# 标题居中
                   panel.grid = element_blank(),  # 移除背景网格线
                   axis.line = element_line(color = "black", size = 0.5)  # 添加x轴和y轴线条
             ) 
)

# 韦恩图
# TCGA差异基因 # log2fc>1,<-1
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/DEG_analyse/step5_差异分析矩阵和结果.Rdata")
# 导出表格
write.csv(expr_for_diff,file = "../出图_24_9_3/附表18_TCGA癌跟癌旁差异分析.csv")
deg_tcga=expr_for_diff[expr_for_diff$change!="nochange",]
tcga=deg_tcga$row
# GSE76250 # allDiff$logFC > 1 & allDiff$adj.P.Val < 0.05
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/geo数据/gse76250_癌跟癌旁分组差异分析_24_9_12.Rdata")
# 导出表格
write.csv(allDiff,file = "../出图_24_9_3/附表19_gse76250差异分析.csv")
deg_gse76250=allDiff[allDiff$type!="not-sig",]
gse76250=deg_gse76250$gene
# 再加上那9个基因取交集
gene=c("UQCRQ","ISCU","FN1","TBX2","NOTCH3","WFDC1","C1S","JUN","PLS3")
# 绘制韦恩图
library(VennDiagram)
# 绘制 3 集合的韦恩图
venn.plot <- venn.diagram(
  x = list(GSE76250 = gse76250, TCGA = tcga, Modelgene = gene),  # 数据
  category.names = c("GSE76250", "TCGA", "Modelgene"),   # 集合标签
  filename = NULL,  # 不保存到文件，直接在R绘图窗口中显示
  col = "transparent",  # 边框颜色
  fill = c("#53A85F", "#E95C59", "#57C3F3"),  # 填充颜色
  # alpha = 0.5,  # 透明度
  cex = 2,  # 字体大小
  cat.cex = 2,  # 标签字体大小
  cat.pos = c(-20, 20, 0),  # 调整标签的位置
  cat.dist = c(0.1, 0.1, 0.1) # 调整标签与图的距离
)

# 显示韦恩图
grid.draw(venn.plot)
dev.off()
#