# 火山图
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/inputgene_24_7_19.Rdata")
# 输入文件
all.markers = markers_genes %>% dplyr::filter(p_val<0.05)
#
# top5= all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10= all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10=as.data.frame(top10)
head(top10)
#   p_val avg_log2FC pct.1 pct.2 p_val_adj cluster    gene
# 1     0   3.413885 0.312 0.023         0    iCAF PLA2G2A
# 2     0   3.086529 0.576 0.208         0    iCAF  CXCL14
# 3     0   3.033794 0.703 0.308         0    iCAF     CFD
# 4     0   2.613782 0.348 0.041         0    iCAF    PI16
# 5     0   2.173342 0.609 0.283         0    iCAF  IGFBP6
# 6     0   2.129693 0.451 0.084         0    iCAF    TNXB
library('ggplot2')
library('dplyr')
library('ggrepel')
library('ggpubr')

#
df=all.markers
# 
head(df)
#         p_val avg_log2FC pct.1 pct.2 p_val_adj cluster    gene
# PLA2G2A     0   3.413885 0.312 0.023         0    iCAF PLA2G2A
# CXCL14      0   3.086529 0.576 0.208         0    iCAF  CXCL14
# CFD         0   3.033794 0.703 0.308         0    iCAF     CFD
# PI16        0   2.613782 0.348 0.041         0    iCAF    PI16
# IGFBP6      0   2.173342 0.609 0.283         0    iCAF  IGFBP6
# TNXB        0   2.129693 0.451 0.084         0    iCAF    TNXB
#添加显著性标签：
df$label <- ifelse(df$p_val_adj<0.05,"adjust P-val<0.05","adjust P-val>=0.05")
head(df)
#
#获取每个cluster中表达差异最显著的10个基因；
table(df$cluster)
top10sigiCAF <- dplyr::filter(df,cluster=="iCAF") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sigiCAF)
#
top10sigmyCAF <- dplyr::filter(df,cluster=="myCAF") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sigmyCAF)
#
top10sigEMT_like <- dplyr::filter(df,cluster=="EMT_like CAF") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sigEMT_like)
#
top10sigVSMC <- dplyr::filter(df,cluster=="VSMC") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sigVSMC)
#
top10sigUndefined_fib <- dplyr::filter(df,cluster=="Undefined_fib") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sigUndefined_fib)
#
top10sigapCAF <- dplyr::filter(df,cluster=="apCAF") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sigapCAF)
#
top10sigUndefined_CAF <- dplyr::filter(df,cluster=="Undefined_CAF") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sigUndefined_CAF)
#
top10sigPericyte <- dplyr::filter(df,cluster=="Pericyte") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sigPericyte)
#
top10sigNAF <- dplyr::filter(df,cluster=="NAF") %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
head(top10sigNAF)

#
#将提取所有cluster的Top10基因表格合并：
top10sig <- rbind(top10sigiCAF,top10sigmyCAF,top10sigEMT_like,top10sigVSMC,top10sigUndefined_fib,
                  top10sigapCAF,top10sigUndefined_CAF,top10sigPericyte,top10sigNAF)
top10sig <- rbind(top10sigmyCAF,top10sigVSMC,top10sigPericyte)
# 筛选df
df <- df %>%
  dplyr::filter(cluster %in% c("myCAF", "VSMC", "Pericyte"))

#新增一列，将Top10的差异基因标记为2，其他的标记为1；
df$size <- case_when(!(df$gene %in% top10sig$gene)~ 1,
                     df$gene %in% top10sig$gene ~ 2)

#提取非Top10的基因表格；
dt <- dplyr::filter(df,size==1)
head(dt)
#
# 然后是绘图第一步！分别绘制需要带geneID标签和不需要带标签的两个火山图并叠加起来。

#绘制每个Cluster Top10以外基因的散点火山图：
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)
p

#叠加每个Cluster Top10基因散点(将散点适当放大强调）：
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p

# 接着开始第二步，画背景的灰色柱子，并把散点叠上去。

#根据图p中log2FC区间确定背景柱长度：
dfbar<-data.frame(x=c("iCAF","myCAF","EMT_like CAF","VSMC","Undefined_fib","apCAF","Undefined_CAF","Pericyte","NAF"),
                  y=c(4,4,4,4,1,5,5,5,4.5))
dfbar<-data.frame(x=c("myCAF","VSMC","Pericyte"),
                  y=c(4,4,5))
# dfbar1<-data.frame(x=c(0,1,2,3,4,5,6,7,8),
#                    y=c(-1.05,-1.1,-1.3,-1.3,-1.8,-1.55,-1.3,-1.9,-0.85))
#绘制背景柱：
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
# +
#   geom_col(data = dfbar1,
#            mapping = aes(x = x,y = y),
#            fill = "#dcdcdc",alpha = 0.6)
p1
#把散点火山图叠加到背景柱上：
pal_npg("nrc",alpha = 0.9)(10)
# [1] "#E64B35E5" "#4DBBD5E5" "#00A087E5" "#3C5488E5" "#F39B7FE5" "#8491B4E5" "#91D1C2E5"
# [8] "#DC0000E5" "#7E6148E5" "#B09C85E5"
p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  # geom_col(data = dfbar1,
  #          mapping = aes(x = x,y = y),
  #          fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)+
  scale_color_manual(values = c("adjust P-val<0.05" = "#E64B35E5", "adjust P-val>=0.05" = "#4DBBD5E5"))
  # scale_color_brewer(palette = "Set2")
  # scale_color_viridis_d(option = "plasma")  # 选择一种 `viridis` 渐变配色
p2
# 最后是第三步！绘制cluster色块并进行叠加：

#添加X轴的cluster色块标签：
# dfcol<-data.frame(x=c(1:9),
#                   y=0,
#                   label=c(0:8))
# mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F","#7E61487F")
#
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
dfcol<-data.frame(x=c(1:3),
                  y=0,
                  label=c(0:2))
mycol <- c("#53A85F","#F3B1A0","#E95C59")

p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.4,
                     color = "black",
                     fill = mycol,
                     alpha = 1, # 调方块颜色深度
                     show.legend = F)
p3

# 现在这张火山图就算初步画完，最后剩下一些geneID标签的添加、主题美化等常规操作。

#给每个Cluster差异表达前Top10基因加上标签：
p4 <- p3+
  geom_text_repel(
    data=top10sig,
    aes(x=cluster,y=avg_log2FC,label=gene),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p4


# 富集
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/inputgene_24_7_19.Rdata")
# 输入文件
all.markers = markers_genes %>% dplyr::filter(p_val<0.05)
#
# top5= all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10= all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10=as.data.frame(top10)
head(top10)
# 筛选下，只要三群
# all.markers_MVP=all.markers %>% filter(cluster==c("myCAF","VSMC","Pericyte"))
table(all.markers$cluster)
# iCAF         myCAF  EMT_like CAF          VSMC Undefined_fib         apCAF 
# 361           291           464           289           116           402 
# Undefined_CAF      Pericyte           NAF 
# 310           188          1099
top150= all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
table(top150$cluster)
all.markers=top150
VSMC=all.markers[all.markers$cluster=='VSMC',]$gene
Pericyte=all.markers[all.markers$cluster=='Pericyte',]$gene
myCAFs=all.markers[all.markers$cluster=='myCAF',]$gene
# total <- list(VSMC=VSMC,Pericyte=Pericyte,myCAFs=myCAFs)

# 多分组富集分析
all.markers_MVP=all.markers %>% dplyr::filter(cluster %in% c("myCAF","VSMC","Pericyte"))
#
library(clusterProfiler)
library(ggplot2)
# 构建分组，转化genesymbol-gene ID。
df_sig=all.markers_MVP

group <- data.frame(gene=df_sig$gene,
                    group=df_sig$cluster)

Gene_ID <- bitr(df_sig$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

#构建文件并分析
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')

# 富集分析，去除冗杂terms！
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
head(data_GO)
dotplot(data_GO, showCategory=10,font.size = 8)
data_Go_fil=data_GO@compareClusterResult

data_GO_sim <- simplify(data_GO,    # 使用 simplify() 函数来简化基因本体（GO）数据，以减少冗余的 GO 术语【相似性大于0.7的term，被认为是冗余的】。
                        cutoff=0.7,
                        by="p.adjust",
                        select_fun=min)


dotplot(data_GO_sim, showCategory=10,font.size = 12)

data_GO_sim_fil <- data_GO_sim@compareClusterResult # 保存条目
save(data_GO_sim_fil,file = "../出图_24_9_3/fig3b三群细胞高表达基因富集.Rdata")

# 导出表格
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/出图_24_9_3/fig3b三群细胞高表达基因富集.Rdata")
write.csv(data_GO_sim_fil,file = "../出图_24_9_3/附表8_三群细胞特征基因富集结果.csv")
# 到ggplot进行可视化编辑
# 每个分组按照p.dajust排序，取top10
top10_per_cluster <- data_GO_sim_fil %>%
  group_by(Cluster) %>%                # 按 Cluster 分组
  slice_min(order_by = p.adjust, n = 10)  # 选取每组 p.adjust 最小的前 10 项
#
df_GO=top10_per_cluster
df_GO$Description <- str_wrap(df_GO$Description, width = 45) # 超过指定长度，自动换行
library(forcats)
df_GO$Description <- as.factor(df_GO$Description)
df_GO$Description <- fct_inorder(df_GO$Description)
#
ggplot(df_GO, aes(Cluster, Description)) +
  geom_point(aes(fill=p.adjust, size=Count), shape=21)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text = element_text(color = 'black', size = 10))+
  scale_fill_gradient(low="red",high="#4DBBD5E5")+
  labs(x=NULL,y=NULL)
ggplot(df_GO, aes(Cluster, Description)) +
  geom_point(aes(color=p.adjust, size=Count), shape=16)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text = element_text(color = 'black', size = 12))+
  scale_color_gradient(low="red",high="blue")+
  labs(x=NULL,y=NULL)
# +
#   coord_flip() # 将图形的 x 轴和 y 轴进行翻转






# wgcna
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/Mime1包/mydata/GSE58812_清洗好的data.Rdata")
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/分型信息.Rdata")
# load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/inputgene_24_7_19.Rdata")
exp_filter=exp
exp_filter[1:4,1:4]

m.mad <- apply(exp_filter,1,mad)
dataExprVar <- exp_filter[which(m.mad >
                                  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
data.mat=t(dataExprVar)
dim(data.mat) # [1]  107 15141
# 要更改为行为样本，列为基因
# 然后只保留MP相关基因
data.mat=data.mat[,colnames(data.mat)%in%geneset1$gene]
dim(data.mat) # [1]  107 319

# 下面过滤异常样本
library(WGCNA)

gsg <- goodSamplesGenes(data.mat, verbose = 3)
gsg$allOK
# TRUE
# 如果返回为true,证明没有缺失值，可以直接进行下一步
# 如果为false,则需要用以下代码进行删除缺失值
# 如果存在异常样本或基因
if (!gsg$allOK) {
  # 异常的基因
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(data.mat)[!gsg$goodGenes], collapse = ",")));
  # 异常的样本
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(data.mat)[!gsg$goodSamples], collapse = ",")));
  # 删除异常的样本和基因
  data.mat = data.mat[gsg$goodSamples, gsg$goodGenes]
}

# 下一步：聚类所有样本，观察是否有离群样本或者异常样本
# 删除剪切线以下的样本，查看图片，选取离群值
# 如果不想删除可以将cutHeight设置高些
# 绘制样本聚类图
sampleTree <- hclust(dist(data.mat), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

# 根据上图，设置一个剪切线，将离群样本删除
# Plot a line to show the cut
abline(h = 25, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 25, minSize = 10)
table(clust)

# 赋值，省的后面改代码
datExpr0=data.mat
outSamples <- datExpr0[(clust==0),]
rownames(outSamples)


# clust 1 contains the samples we want to keep.
keepSamples <- (clust==1)
datExpr <- datExpr0[keepSamples, ]
sampleTree2 <- hclust(dist(datExpr), method = "average");
par(cex = 0.5);
par(mar = c(0,6,2,0))
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


# 加入临床特征，前面准备好了
# 有一点需要注意，这里能作为输入的必须为数值型特征。
# 把high改为1，low改为0
head(a)
# a=a[,-1]
#            MP_score MP_group
# GSM1419942 2.548510     high
# GSM1419943 2.216109      low
# GSM1419944 2.081702      low
# GSM1419945 2.649311     high
# GSM1419946 2.475124     high
# GSM1419947 2.511488     high
a$group=NA
a$group[which(str_detect(a$MP_group, "^high"))] <- "1"
a$group[which(str_detect(a$MP_group, "^low"))] <- "0"
table(a$group)
# 0  1 
# 54 53
datTraits=a
datTraits[,3]=as.numeric(datTraits[,3])
head(datTraits)
#            MP_score MP_group group
# GSM1419942 2.548510     high     1
# GSM1419943 2.216109      low     0
# GSM1419944 2.081702      low     0
# GSM1419945 2.649311     high     1
# GSM1419946 2.475124     high     1
# GSM1419947 2.511488     high     1
#
identical(rownames(datExpr),rownames(datTraits))
#
datTraits=datTraits[rownames(datTraits)%in%rownames(datExpr),]
identical(rownames(datExpr),rownames(datTraits))
# [1] FALSE
# datExpr=datExpr[rownames(datExpr)%in%rownames(datTraits),]
# identical(rownames(datExpr),rownames(datTraits))
# # [1] TRUE
datTraits=datTraits[3]
rownames(datTraits)
# collectGarbage()
# 增加性状信息后，再次聚类
sampleTree2=hclust(dist(datExpr),method = "average")
# 绘图
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr0,datExpr, datTraits, file = "../../Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/GEO-01-dataInput.RData")


# 2、网路构建和模块检测
# 2.1 选择软阈值
# 使用函数**pickSoftThreshold()**选择适当的软阈值。
rm(list = ls())
options(stringsAsFactors = F)

library(WGCNA)

load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/GEO-01-dataInput.RData")

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
par(mar = c(3,6,2,1))
par(mfrow = c(1,2))
cex1 = 0.9 #可以修改
# 拟合指数与power值散点图，无标度拓扑拟合指数
# 一般选择在0.9以上的，第一个达到0.9以上的数值

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", # n表示不绘图
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")

# best
sft$powerEstimate
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red") # 这边h设置为cex1的值

# 越平滑越好
# 平均连通性与power值散点图
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=10,col="red")
# 其中sft$powerEstimate是推荐的最优软阈值，最优值为4。左图y轴是无标度拓扑拟合指数，右图y轴是平均连通度。
# 其中横坐标为软阈值的梯度，第一幅图的纵坐标为无标度网络适应系数，越大越好；第二幅图的纵坐标为节点的平均连通度，越小越好。
sft$powerEstimate
## 如果推荐的最优软阈值为NA，则表面系统无法给出合适的软阈值，这时候就需要自己挑选软阈值。
sft
softpower=sft$powerEstimate
adjacency=adjacency(datExpr,power = softpower)
softpower


# TOM矩阵
TOM=TOMsimilarity(adjacency)
dissTOM=1-TOM
# 基因聚类
geneTree=hclust(as.dist(dissTOM),method = "average")
plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity",
     labels=FALSE,hang=0.04)
# 动态剪切模块识别
minModuleSize=20 # 最小单个模块包含的基因数
dynamicMods=cutreeDynamic(dendro = geneTree,distM = dissTOM,
                          deepSplit = 2,pamRespectsDendro = FALSE,
                          minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors=labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree,dynamicColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE,hang=0.1,
                    addGuide = TRUE,guideHang = 0.1,
                    main="Gene dendrogram and module colors")

net <- blockwiseModules(datExpr, power = 3,
                        maxBlockSize = 2000,# 最大模块数量
                        TOMType = "unsigned", 
                        minModuleSize = 20, #用于模块检测的最小模块尺
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.25, # 用于模块合并的树形图切割高度
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "../../../Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/TNBCTOM",
                        verbose = 3)
# minModuleSize表示用于模块检测的最小模块尺寸。mergeCutHeight表示用于模块合并的树形图切割高度。这两个值越大，模块越少。saveTOMFileBase用来设置数据保存位置及名称。

table(net$colors)
# 0   1 
# 91 228 
# 查看划分的模块数和每个模块里面包含的基因个数。0代表无法识别的基因数。
# 用于模块识别的分层聚类树形图存储在net$dendprograms[[1]]中。可以与模型颜色分配一起显示。
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.1,
                    addGuide = TRUE, guideHang = 0.04)
# 除此之外，还有分布网络构建和模块检测，以及处理大数据集的：分块网络构建和模块检测。

# 保存参数
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs;
geneTree <- net$dendrograms[[1]];

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "../../Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/GEO-02-networkConstruction-auto.RData")




# 3、将模块与表型数据关联并识别重要基因
# 3.1 模块与性状关系图
rm(list = ls())
options(stringsAsFactors = F)

library(WGCNA)

load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/GEO-01-dataInput.RData")
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/GEO-02-networkConstruction-auto.RData")

# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 10, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
# 图中展示了3个模块与1种性状的关系，其中红色代表模块与性状正相关，蓝色则代表负相关。
# 我们要选择与性状最相关的模块，那么就是MEturquoise对应的模块，其P值为0.65。

# 3、青色模块与MP评分的关系（聚类图和热图）
# 使用plotEigengeneNetworks函数再次验证青色模块与肿瘤是否有关联。
# Recalculate module eigengenes
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
# 注意：记得修改代码，选择你感兴趣的性状
samples <- as.data.frame(datTraits$group);
names(samples) <- "samples"
# Add the weight to existing module eigengenes
MET <- orderMEs(cbind(MEs, samples))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,9)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(4,5,1,2),
                      cex.lab = 0.8, xLabelsAngle = 90)
# 聚类图中可以看出samples与MEturquiose高度相关，在热图中青色模块与samples也有显著相关性。说明青色模块就是与肿瘤最相关的模块。
# 也可以将聚类图和热图分开展示。
# Plot the dendrogram
sizeGrWindow(6,6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)



# 4、青色模块中MP的GS与MM的相关性（散点图）
# 我们将基因显著性 (Gene Significance, GS) 定义为基因与性状之间相关性的（绝对值），以此量化单个基因与我们感兴趣的性状（权重）之间的关联。
# 对于每个模块，我们还定义了一个定量指标模块成员(module membership, MM），即模块特征基因与基因表达谱的相关性。这样我们就可以量化阵列上所有基因与每个模块的相似性。代码实现如下：
samples <- as.data.frame(datTraits$group);
names(samples) <- "samples"
# names (colors) of the modules
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")
geneTraitSignificance <- as.data.frame(cor(datExpr, samples, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(samples), sep="")
names(GSPvalue) <- paste("p.GS.", names(samples), sep="")
# 注意：其中前两行用来选择你感兴趣的性状，记得修改成对应的数据。
# 针对感兴趣的模块，将GS与MM值可视化，代码如下：
modNames
# [1] "blue"      "grey"      "turquoise"
module <- "blue"
column <- match(module, modNames);
moduleGenes <- moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
# [1] "#E64B35E5" "#4DBBD5E5" "#00A087E5" "#3C5488E5" "#F39B7FE5" "#8491B4E5" "#91D1C2E5"
# [8] "#DC0000E5" "#7E6148E5" "#B09C85E5"
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for sample type",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "#029AE5E5")
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for sample type",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "turquoise")

# 注意：module对应的是你感兴趣模块的颜色，记得修改图中的文字。
# 再注意：verboseScatterplot函数中的col本来是“col = module”，但是由于黄色太不明显了，我就用另一种黄色代替了。

# 图中的每一个点代表一个基因。显然，图中GS 和 MM 高度相关，说明与肿瘤高度显著相关的基因往往也是黄色模块中最重要（核心）的元素。

# 
# 5、保存与肿瘤最相关的基因
colnames(datExpr)
table(moduleGenes)
sel_genes1 <- colnames(datExpr)[moduleColors=="turquoise"]
sel_genes2 <- colnames(datExpr)[moduleColors=="blue"]
# write.table(sel_genes, file = "../../../../Breastcancer_EMBOJ/emboj_预后模型/wgcna/module_gene_name.txt", quote = F, 
#             row.names = F,col.names = F)
sel_genes=c(sel_genes1,sel_genes2)
save(sel_genes,sel_genes1,sel_genes2,file = "../../Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/最相关模块基因.Rdata")





# 单因素cox
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/wgcna/new/最相关模块基因.Rdata")
# 单因素回归
# 通过一个for循环对所有目标基因进行回归分析，并且以dataframe的形式对结果进行输出：
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/Mime1包/mydata/GSE58812_不分割log的清洗好的data.Rdata")
Train[1:4,1:4]
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
# 列名不能自行修改，如果完全不懂代码的小白，最好完全照着我的数据，保持列名与我的一致，不然容易报错。
rt <- data
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


# 设置表格内容
tabletext <- cbind(
  c("Gene", outTab$id),
  c("HR", sprintf("%.3f", outTab$HR)),
  c("95% CI", paste(sprintf("%.3f", outTab$HR.95L), "-", sprintf("%.3f", outTab$HR.95H))),
  c("P value", sprintf("%.3e", outTab$pvalue))
)
forestplot(
  tabletext,
  mean = c(NA, outTab$HR), 
  lower = c(NA, outTab$HR.95L), 
  upper = c(NA, outTab$HR.95H),
  zero = 1,  # 指定零点线
  boxsize = 0.3,  # 方块大小
  lineheight = unit(1.5, "cm"),  # 每一行的高度
  col = fpColors(box = color_palette, lines = "#6666FF", zero = "gray50"), # 设置颜色：显著性HR点用 `color_palette` 动态分配，线条为深蓝色，零点线为灰色
  xlog = FALSE,  # 对数尺度展示
  ci.vertices = TRUE,  # 显示置信区间顶点
  ci.vertices.height = 0.1,  # 置信区间顶点高度
  title = "Enhanced Forest Plot of HR with 95% CI",
  graphwidth = unit(8, "cm"),  # 增加图形宽度
  xticks = c(0,0.5, 1, 2, 4, 8,10,12),  # 自定义X轴刻度
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 12, fontface = "bold"),  # 基因名称加粗
    ticks = gpar(fontsize = 10),  # 刻度字体大小
    xlab = gpar(fontsize = 12)  # X轴标签字体
  ),
  hrzl_lines = gpar(col = "#444444"),  # 水平线颜色
  grid = gpar(lty = 2, col = "#CCCCCC"),  # 添加背景网格线
  clip = c(0, 10),  # 设置置信区间的上下限
  is.summary = FALSE,  # 取消汇总线
  new_page = TRUE  # 在新页面上绘制
)

forestplot(
  tabletext,  # 包含文本的矩阵，用于显示森林图左侧的基因信息、HR值、95% CI 和 p 值等
  mean = c(NA, outTab$HR),  # 中心点，显示HR值，第一项NA是为了对齐表头
  lower = c(NA, outTab$HR.95L),  # 置信区间的下限值
  upper = c(NA, outTab$HR.95H),  # 置信区间的上限值
  zero = 1,  # 指定零点线为1（HR=1是分界点，表示无效假设）
  boxsize = 0.3,  # 设置方块大小，表示HR的可视化点
  lineheight = unit(1.5, "cm"),  # 每一行的高度，确保表格和点之间有足够的间距
  col = fpColors(box = color_palette, lines = "darkblue", zero = "gray50"),  # 设置颜色：显著性HR点用 `color_palette` 动态分配，线条为深蓝色，零点线为灰色
  xlog = FALSE,  # 设置X轴为线性尺度而不是对数尺度
  ci.vertices = TRUE,  # 在置信区间末端显示顶点
  ci.vertices.height = 0.1,  # 设置置信区间顶点的高度
  title = "Enhanced Forest Plot of HR with 95% CI",  # 设置森林图的标题
  graphwidth = unit(8, "cm"),  # 调整图形的宽度为8厘米
  xticks = c(0,0.5, 1, 2, 4, 8,10,11),  # 设置X轴刻度范围，范围从0到5，步长为0.5
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 12, fontface = "bold"),  # 基因名称字体大小为12，并加粗
    ticks = gpar(fontsize = 10),  # X轴刻度字体大小
    xlab = gpar(fontsize = 12)  # X轴标签字体大小
  ),
  hrzl_lines = gpar(col = "#444444"),  # 水平线的颜色设置为灰色
  grid = gpar(lty = 2, col = "#CCCCCC"),  # 添加背景网格线，颜色为浅灰色，虚线样式
  clip = c(0, 5),  # 限制置信区间的显示范围，最小为0，最大为5
  is.summary = FALSE,  # 指定该图不是汇总图
  new_page = TRUE  # 在新页面上开始绘图
)

