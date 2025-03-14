distribution_OR <- function(
    meta_data,
    celltype_column,
    celltype_level = NULL,
    condition_column,
    condition_level = NULL
){
  library(tidyverse)

  colnames(meta_data)[which(colnames(meta_data) == celltype_column)] = "celltypE"
  colnames(meta_data)[which(colnames(meta_data) == condition_column)] = "conditioN"
  
  if(is.null(celltype_level)){
    meta_data$celltypE = as.character(meta_data$celltypE)
    meta_data$celltypE = factor(meta_data$celltypE,levels = sort(unique(meta_data$celltypE)))
  } else {
    meta_data$celltypE = factor(meta_data$celltypE,levels = celltype_level)
  }
  
  if(is.null(condition_level)) {
    meta_data$conditioN = as.character(meta_data$conditioN)
    meta_data$conditioN = factor(meta_data$conditioN,levels = sort(unique(meta_data$conditioN)))
  } else {
    meta_data$conditioN = factor(meta_data$conditioN,levels = condition_level)
  }

  count.dist = as.data.frame(table(meta_data$celltypE,meta_data$conditioN))
  count.dist = spread(count.dist,key = Var2,value = Freq)
  rownames(count.dist) = count.dist$Var1
  count.dist$Var1 =NULL
  count.dist = as.matrix(count.dist)

  ######################
  #
  library(data.table)
  
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.DT <- as.data.frame(count.dist)
  setDT(count.dist.DT,keep.rownames=T)
  
  count.dist.DT.melt <- data.table::melt(count.dist.DT,id.vars="rn")
  colnames(count.dist.DT.melt) <- c("rid","cid","count")
  
  
  library(plyr)
  test.res <- as.data.table(
    ldply(
      seq_len(nrow(count.dist.DT.melt)), function(i){
        this.row <- count.dist.DT.melt$rid[i]
        this.col <- count.dist.DT.melt$cid[i]
        this.c <- count.dist.DT.melt$count[i]
        
        this.m <- matrix(
          c(this.c,
            sum.row[this.row]-this.c,
            sum.col[this.col]-this.c,
            sum(sum.col)-sum.row[this.row]-sum.col[this.col]+this.c),
          ncol=2)
        #        #                  this celltype | not this celltype
        #     this tissue|      a       |         c
        #not this tissue|      b       |         d
        
        
        #阈值不固定
        tmp.res <- fisher.test(this.m)
        data.frame(
          rid=this.row,
          cid=this.col,
          p.value=tmp.res$p.value,
          OR=tmp.res$estimate # 约为 a*d / b*c
        )
      }
    )
  )
  
  test.res <- merge(count.dist.DT.melt,test.res,by=c("rid","cid"))
  test.res[,adj.p.value:=p.adjust(p.value,"BH")]
  test.res = as.data.frame(test.res)

  #############################################################
  dist.p <- reshape2::dcast(test.res,rid~cid,value.var="p.value")
  dist.OR <- reshape2::dcast(test.res,rid~cid,value.var="OR")
  dist.p.adj <- reshape2::dcast(test.res,rid~cid,value.var="adj.p.value")
  rownames(dist.p) = dist.p$rid
  rownames(dist.OR) = dist.OR$rid
  rownames(dist.p.adj) = dist.p.adj$rid
  dist.p$rid=NULL
  dist.OR$rid=NULL
  dist.p.adj$rid=NULL
  
  return(list(
      "dist.p"=dist.p,
      "dist.OR"=dist.OR,
      "dist.p.adj"=dist.p.adj))
}








#### cellchat
sce=CreateSeuratObject(counts = sce@assays$RNA@counts,meta.data = sce@meta.data)
DefaultAssay(sce)
#
scRNAlist=sce
DimPlot(scRNAlist,reduction = "tsne",label = T)
table(scRNAlist@active.ident)
table(scRNAlist@meta.data$forth_integrated_bcell_6_17)  ## 28种细胞类型
Idents(scRNAlist)=scRNAlist$forth_integrated_bcell_6_17

data.input  <- scRNAlist@assays$RNA@data

table(scRNAlist@meta.data$Tissue)
# ER   HER2 Normal   TNBC 
# 67668  46078  54614  21951 

meta = scRNAlist@meta.data # a dataframe with rownames containing cell mata data
unique(meta$forth_integrated_bcell_6_17) # check the cell labels 也是28种细胞类型
# cell_type就是labels


colnames(meta)[17]="labels"

## 创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

## 设置配体受体交互数据库
## 我们的数据库 CellChatDB 是一个手动整理的文献支持的配体受体在人和小鼠中的交互数据库。
## 小鼠中的CellChatDB包含2，021个经验证的分子相互作用，包括60%的自分泌/旁分泌信号相互作用、21%的细胞外基质（ECM）受体相互作用和19%的细胞-细胞接触相互作用
## 人的CellChatDB包含1，939个经验证的分子相互作用，包括61.8%的自分泌/旁分泌信号相互作用、21.7%的细胞外基质（ECM）受体相互作用和16.5%的细胞-细胞接触相互作用
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
colnames(CellChatDB$interaction)
table(CellChatDB$interaction$annotation)
# Cell-Cell Contact       ECM-Receptor Secreted Signaling 
# 319                421               1199 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
unique(CellChatDB$interaction$annotation)
# [1] "Secreted Signaling" "ECM-Receptor"      
# [3] "Cell-Cell Contact" 
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use  ##




##预处理用于细胞通信分析的表达数据
cellchat <- subsetData(cellchat) # 取出表达数据
# cellchat <- subsetData(cellchat，features = NULL) # 这边可以选择感兴趣的基因
## 不想改成多线程
future::plan("multisession", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)  ## 寻找高表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)  ## 寻找高表达的通路
future::plan("multisession", workers = 1)
# cellchat <- projectData(cellchat, PPI.human) ##投影倒PPI
cellchat <- projectData(cellchat, PPI.human)
## 我们这边是基因，但是有的时候基因跟蛋白质不一定对得上，还有一个可能是基因有的不会翻译成蛋白质
## 这边做了一个投影到PPI相当于做了一个翻译的过程

# 第二部分
## 细胞通信网络的推断
## 计算通信概率并推断cellchat网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
?computeCommunProb
## 这边raw.use 如果选FALSE,则这一步是用上面计算的蛋白质的结果计算
## 如果raw.use选TRUE，则这一步是用基因的矩阵计算
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
## 在computeCommunProb中，我们提供了一个选项，用于使用其他方法
## 默认情况下，CellChat 使用一种统计学上强大的均值方法，称为"trimean"，
## "trimean"大约是25%的截断平均值，这意味着如果一组表达细胞的百分比低于25%，则平均基因表达为零【为0的这些基因就不会再纳入细胞通讯的一个计算】。
## 要使用 10% 截断的平均值，用户可以设置type = "truncatedMean"和对trim = 0.1。

cellchat <- filterCommunication(cellchat, min.cells = 10) #去掉通讯数量很少的细胞
df.net <- subsetCommunication(cellchat)  ##将细胞通讯预测结果以数据框的形式取出
# df.net <- subsetCommunication(cellchat)  ##返回一个数据框架，该数据框架由配体/受体级别的所有推断细胞通信组成。设置slot.name = "netP"可以在信号通路级别访问推断的通信
# df.netp <- subsetCommunication(cellchat,solt.name="netP") ##只取通路，数据结构更简单
levels(cellchat@idents)


# 推测的每个"配体-受体"对的细胞间通信网络和每个"信号通路"分别存储在“net”和“netP”槽中。


## 在信号通路级别推断细胞-细胞通信
## CellChat 通过总结与每个信号通路相关的所有配体-受体相互作用的通信概率，来计算信号通路级别上的通信概率。
## NB：每个配体受体对和每个信号通路的推断细胞间通信网络分别存储在插槽"net"和"netP"中。
cellchat <- computeCommunProbPathway(cellchat)
# 计算整合的细胞通信网络
# 我们可以通过计算链接数或汇总通信概率来计算整合的细胞通信网络。用户还可以通过设置sources.use和targets.use`
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways   ##看看有哪些信号通路
# [1] "COLLAGEN"  "LAMININ"   "FN1"       "MIF"       "CD99"      "MK"        "CXCL"     
# [8] "APP"       "THBS"      "ANGPTL"    "CCL"       "NOTCH"     "PTN"       "TENASCIN" 
# [15] "FGF"       "MPZ"       "GAS"       "ESAM"      "PERIOSTIN" "SEMA3"     "GRN"      
# [22] "CD46"      "PDGF"      "PROS"      "GALECTIN"  "ADGRE5"    "ANGPT"     "SELL"     
# [29] "NT"        "TWEAK"     "TGFb"      "CADM"      "SEMA5"     "JAM"       "ITGB2"    
# [36] "EGF"       "HSPG"      "BMP"       "NEGR"      "IGF"       "AGRN"
head(cellchat@LR$LRsig)  ##看看具体的配受体情况
# 我们还可以可视化整合的细胞通信网络。例如，使用圆图显示任意两个细胞组之间的相互作用次数或总交互强度（比重）。
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)  ##设置图片的一个布局
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")  ##互作的数量
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Interaction weights/strength")  ##互作的强度

#

## 互作数量和重要性图
## 由于细胞通信网络复杂，
## 我们可以检查每个细胞组发送的信号。在这里，我们还控制参数edge.weight.max，以便我们可以比较不同网络之间的边缘权重。
#
# load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/cellchat/cellchat_all/cellchat_all.Rdata") # 这是原来拿所有组计算的，包括er,her2,tnbc,normal


# KM
#
rm(list=ls())
gc()
# 
# 
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/ssgsea/new/cellMarker_ssGSEA_全部乳腺癌亚型_24_7_19.Rdata") # 制备的cellMarker
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/Mime1包/mydata/GSE58812_清洗好的data.Rdata")
expr=as.matrix(exp)
gsva_data <- gsva(expr,cellMarker, method = "ssgsea",abs.ranking = TRUE)  # 默认参数kcdf = "Gaussian"，适用于对数转换的microarray、RNA-seq的log-CPMs、log-RPKMs或log-TPMs。当输入表达的矩阵是RNA-seq的raw Count时，这个参数应该设置为kcdf = "Poisson"。

a <- gsva_data %>% t() %>% as.data.frame()
results=a
# 用我们整理好的预后信息
head(cli)
# 使用sub函数提取冒号后面的值
cli$OS <- sub(".*: ", "", cli$characteristics_ch1.4)
cli$OS.time <- sub(".*: ", "", cli$characteristics_ch1.5)
#
cli=cli[,c(3,4)]
head(cli)
#
# sub 函数用于字符串替换。
# "*.: " 匹配任何字符直到最后一个冒号和空格。
# "" 替换匹配的部分为空字符串，即删除匹配的部分。
sur_data=cli
sur_data$OS.time=as.numeric(sur_data$OS.time)
sur_data$OS=as.numeric(sur_data$OS)
# #去掉生存信息不全或者生存时间小于30天的样本，样本纳排标准不唯一，且差别很大.
sur_data =sur_data[sur_data$OS.time >= 30,]
sur_data = sur_data[!(is.na(sur_data$OS.time)|is.na(sur_data$OS)),]
sur_data[1:4,1:2]

results[1:4,1:4] 
results2=as.data.frame(results)

sur_data2=sur_data[rownames(sur_data) %in% rownames(results2),]
results3=results2[rownames(results2) %in% rownames(sur_data2),]
results3$bcr_patient_barcode=rownames(results3)
sur_data2$bcr_patient_barcode=rownames(sur_data2)
#
aimplot <- left_join(x=sur_data2,y=results3,by='bcr_patient_barcode')
save(aimplot,file = "../../../Breastcancer_EMBOJ/emboj_预后模型/ssgsea/new/aimplot_细胞评分结果.Rdata")

# 最佳分割点
colnames(aimplot)[7]="EMT_likeCAF"
res.cut <- surv_cutpoint(aimplot, time = "OS.time", event = "OS",
                         variables = c("Pericyte",      
                                       "myCAF","VSMC","CD8_STMN1",
                                       "Macrophage","iCAF",             
                                       "EMT_likeCAF","Edothelials","cDC3" ,              
                                       "Plasma_Bcells","apCAF",              
                                       "Myepithelials","cDC2",               
                                       "cDC1","Memory_Bcells","CD4_Tem",            
                                       "CD4_Treg","Monocyte","Mast",               
                                       "CD8_Teff","CD4_CXCL13","NK",                 
                                       "NAF","Tact_IFI6","Naive_Bcells",       
                                       "CD4_HSPA1A","pDC"   )) # 
summary(res.cut)
# 使用surv_categorize()函数根据最佳截断值对数据进行分组（高表达/低表达组）：
res.cut2 <- surv_categorize(res.cut)
head(res.cut2)
table(res.cut2$Pericyte)
table(res.cut2$myCAF)
table(res.cut2$VSMC)
table(res.cut2$CD8_STMN1)

# 配色
# RColorBrewer中的调色板
palette = "Set2"  # Set1, Set2, Paired 等都是不错的选择
# 使用 ggsci 的科学期刊配色：
library(ggsci)
palette = "npg"  # 使用 Nature 期刊风格配色
# 根据最佳截断点绘制生存曲线：
fit <- survfit(Surv(OS.time, OS) ~ Pericyte, data = res.cut2)
ggsurvplot(fit, data = res.cut2,
           conf.int = TRUE,  pval = TRUE,
           surv.median.line = "hv",
           risk.table = TRUE, palette = "hue")
ggsurvplot(fit, data = res.cut2,
           conf.int = TRUE,  pval = TRUE,
           surv.median.line = "hv", # 中位生存线
           # risk.table = TRUE, 
           palette = c("blue","purple")
           
)
ggsurvplot(fit, data = res.cut2,
           conf.int = TRUE,  
           pval = TRUE,
           surv.median.line = "hv",
           # risk.table = TRUE, 
           # risk.table.height = 0.25,  # 调整风险表高度
           # risk.table.col = "strata", # 风险表按分层显示颜色
           palette = "npg",           # 配色方案
           ggtheme = theme_minimal() # 使用干净的主题
           )  
#


