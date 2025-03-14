gene <- dd$gene## 转换
library(clusterProfiler)
gene = bitr(gene, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db')## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=dd$logFC,
                      SYMBOL = dd$gene)
gene_df <- merge(gene_df,gene,by='SYMBOL')
head(gene_df)
#     SYMBOL       logFC ENTREZID
# 1     A1BG  0.06575477        1
# 2 A1BG-AS1 -0.02881001   503538
# 3     A1CF -0.05084698    29974
# 4      A2M  0.06471498        2
# 5  A2M-AS1  0.13557040   144571
# 6    A2ML1  0.10653129   144568
## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)


# 5.运行GSEA分析
# 
# 从GESA(https://www.gsea-msigdb.org/gsea/downloads.jsp)的官网上，下载一个gmt文件

library(clusterProfiler)
## 读入hallmarks gene set，从哪来？ 这边要下载entrezid版本
hallmarks <- read.gmt('../../new_analyse_24_4_17/gsea/h.all.v2023.2.Hs.entrez.gmt')
# 需要网络
y <- GSEA(geneList,TERM2GENE =hallmarks)
# 作图看整体分布

### 看整体分布library(ggplot2)
dotplot(y,showCategory=12,split='.sign')+facet_grid(~.sign)
dotplot(y,showCategory=50,split='.sign')+facet_grid(~.sign)+
  
  # 调整字体大小
  theme(
    axis.title = element_text(size = 14),  # 坐标轴标题字体大小
    axis.text = element_text(size = 12),   # 坐标轴刻度标签字体大小
    strip.text = element_text(size = 14),  # 分面标签字体大小
    legend.title = element_text(size = 12),  # 图例标题字体大小
    legend.text = element_text(size = 10)    # 图例文本字体大小
  )
# 6.特定通路作图

yd <- data.frame(y)

library(enrichplot)
gseaplot2(y,'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',color = 'red',pvalue_table = T)
gseaplot2(y,'HALLMARK_TGF_BETA_SIGNALING',color = 'red',pvalue_table = T)
gseaplot2(y,'HALLMARK_NOTCH_SIGNALING',color = 'red',pvalue_table = T)
gseaplot2(y,'HALLMARK_WNT_BETA_CATENIN_SIGNALING',color = 'red',pvalue_table = T)
gseaplot2(y,'HALLMARK_PI3K_AKT_MTOR_SIGNALING',color = 'red',pvalue_table = T)
gseaplot2(y,'HALLMARK_GLYCOLYSIS',color = 'red',pvalue_table = T)
gseaplot2(y,'HALLMARK_OXIDATIVE_PHOSPHORYLATION',color = 'red',pvalue_table = T)
gseaplot2(y,'HALLMARK_OXIDATIVE_PHOSPHORYLATION',color = 'red',pvalue_table = T)

##对于多个通路绘制在一起：：：
pathway=c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING","HALLMARK_NOTCH_SIGNALING","HALLMARK_GLYCOLYSIS","HALLMARK_OXIDATIVE_PHOSPHORYLATION")
gseaplot2(y,pathway,color = c('red',"blue"),pvalue_table = F)

#
save(yd,y,file = "../101机器学习/得到模型后续分析/GSEA_riskscore_用的是gse58812去计算_24_8_30.Rdata")











# hallmark相关性
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/不用封装的跑/选定模型计算好的riskscore_24_8_29.Rdata")
#
# 想这么做，就需要先通过gsva打分的方法计算每个样本各自通路的评分
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/得到模型后续分析/gse58812_重新清洗cli_24_8_29.Rdata")
exp_filter=exp
exp_filter[1:4,1:4]
#               TCGA-BH-A0E0-11A TCGA-BH-A18V-11A TCGA-BH-A1FC-11A TCGA-E2-A158-11A
# RP11-368I23.2        0.5770281       0.00000000       0.33007291       0.93334894
# RP11-742D12.2        0.0000000       0.03483256       0.02550387       0.02258519
# RAB4B                3.0713037       3.90911352       3.50316793       2.76627551
# AC104183.2           0.0000000       0.00000000       0.00000000       0.00000000
expr=as.matrix(exp_filter)

# 2，获取目标基因集
# 根据自己的需要选择MSigDB数据库中的基因集
# 2.1 手动下载
# 
# 进入http://www.gsea-msigdb.org/gsea/msigdb/index.jsp后选择需要下载的基因集，然后使用R读取下载好的gmt格式的文件。
# 下载50个肿瘤特征基因集合
# 多取几个基因集，取交集
# data1=clusterProfiler::read.gmt("../101机器学习/得到模型后续分析/h.all.v2024.1.Hs.symbols.gmt") # 返回的数据框
# geneset=GSA::GSA.read.gmt("new_analyse_24_4_17/h.all.v2023.2.Hs.symbols.gmt")
# geneset=mogsa::prepMsigDB("new_analyse_24_4_17/h.all.v2023.2.Hs.symbols.gmt") #list
geneset <- cogena::gmt2list("../101机器学习/得到模型后续分析/h.all.v2024.1.Hs.symbols.gmt") #list
# $HALLMARK_PROTEIN_SECRETION
# [1] "ABCA1"    "ADAM10"   "ANP32E"   "AP1G1"    "AP2B1"    "AP2M1"    "AP2S1"    "AP3B1"   
# [9] "AP3S1"    "ARCN1"    "ARF1"     "ARFGAP3"  "ARFGEF1"  "ARFGEF2"  "ARFIP1"   "ATP1A1"  
# [17] "ATP6V1B1" "ATP6V1H"  "ATP7A"    "BET1"     "BNIP3"    "CAV2"     "CD63"     "CLCN3"   
# [25] "CLN5"     "CLTA"     "CLTC"     "COG2"     "COPB1"    "COPB2"    "COPE"     "CTSC"    
# [33] "DNM1L"    "DOP1A"    "DST"      "EGFR"     "ERGIC3"   "GALC"     "GBF1"     "GLA"     
# [41] "GNAS"     "GOLGA4"   "GOSR2"    "ICA1"     "IGF2R"    "KIF1B"    "KRT18"    "LAMP2"   
# [49] "LMAN1"    "M6PR"     "MAPK1"    "MON2"     "NAPA"     "NAPG"     "OCRL"     "PAM"     
# [57] "PPT1"     "RAB14"    "RAB22A"   "RAB2A"    "RAB5A"    "RAB9A"    "RER1"     "RPS6KA3" 
# [65] "SCAMP1"   "SCAMP3"   "SCRN1"    "SEC22B"   "SEC24D"   "SEC31A"   "SGMS1"    "SH3GL2"  
# [73] "SNAP23"   "SNX2"     "SOD1"     "SSPN"     "STAM"     "STX12"    "STX16"    "STX7"    
# [81] "TMED10"   "TMED2"    "TMX1"     "TOM1L1"   "TPD52"    "TSG101"   "TSPAN8"   "USO1"    
# [89] "VAMP3"    "VAMP4"    "VAMP7"    "VPS45"    "VPS4B"    "YIPF6"    "YKT6"     "ZW10"    
# 
# $HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
# [1] "ABCC1"   "ATOX1"   "CAT"     "CDKN2D"  "EGLN2"   "ERCC2"   "FES"     "FTL"     "G6PD"   
# [10] "GCLC"    "GCLM"    "GLRX"    "GLRX2"   "GPX3"    "GPX4"    "GSR"     "HHEX"    "HMOX2"  
# [19] "IPCEF1"  "JUNB"    "LAMTOR5" "LSP1"    "MBP"     "MGST1"   "MPO"     "MSRA"    "NDUFA6" 
# [28] "NDUFB4"  "NDUFS2"  "NQO1"    "OXSR1"   "PDLIM1"  "PFKP"    "PRDX1"   "PRDX2"   "PRDX4"  
# [37] "PRDX6"   "PRNP"    "PTPA"    "SBNO2"   "SCAF4"   "SELENOS" "SOD1"    "SOD2"    "SRXN1"  
# [46] "STK25"   "TXN"     "TXNRD1"  "TXNRD2" 

# 
# 2.2 msigdbr包
# 
# 直接使用msigdbr包内置好的基因集，含有多个物种 以及 多个基因集，通过参数选择物种以及数据集，较为方便。推荐！

# library(msigdbr)
# msigdbr_species() #列出有的物种

#选择基因集合
# ?msigdbr
# human_KEGG = msigdbr(species = "Homo sapiens", #物种
#                      category = "C2",
#                      subcategory = "KEGG") %>% 
#   dplyr::select(gs_name,gene_symbol)#这里可以选择gene symbol或者ID
# human_KEGG_Set = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)#list
#

# A：如果你的研究是其中的物种就可以无缝做GSEA 和 GSVA了。
# 
# B：如果研究的物种不在其中，也可以自定义基因集，注意转为对应的形式。human_KEGG_Set 为基因集合的列表形式。

# 二 GSVA分析
# 
# 1, GSVA分析
# 
# 数据准备好后，加载GSVA包，一个gsva函数就可以得到GSVA的结果了。

library(GSVA)
gsva.kegg <- gsva(expr, gset.idx.list = geneset, 
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=1)
head(gsva.kegg)
#                                GSM1419942  GSM1419943  GSM1419944  GSM1419945
# HALLMARK_ADIPOGENESIS         0.319234180 -0.02995030 -0.12036696  0.12704217
# HALLMARK_ALLOGRAFT_REJECTION -0.373894018  0.09860410 -0.45872752  0.29896703
# HALLMARK_ANDROGEN_RESPONSE    0.375729281  0.13031497 -0.02026651  0.14336347
# HALLMARK_ANGIOGENESIS         0.073590349  0.15138661 -0.23981522  0.44516777
# HALLMARK_APICAL_JUNCTION     -0.006224164 -0.13768544 -0.15815801  0.19611173
# HALLMARK_APICAL_SURFACE      -0.100360665 -0.08079993 -0.29144559 -0.03340218
# 
# 行为目标基因集，列为celltype ，数值为gsva分数。
# 
# 这里需要注意，如果输入矩阵为log转化后的连续表达矩阵指则设置kcdf参数为"Gaussian"，如果是counts矩阵则设置kcdf为"Poisson"。
save(gsva.kegg,file = "../101机器学习/得到模型后续分析/gse数据_hallmark通路富集评分_24_8_30.Rdata")
# 2, 绘制热图
# 以结果的前50个绘制示例热图，可以自选择重点的通路
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/不用封装的跑/选定模型计算好的riskscore_24_8_29.Rdata")
#
# 想这么做，就需要先通过gsva打分的方法计算每个样本各自通路的评分
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/得到模型后续分析/gse58812_重新清洗cli_24_8_29.Rdata")
exp_filter=exp
exp_filter[1:4,1:4]
#               TCGA-BH-A0E0-11A TCGA-BH-A18V-11A TCGA-BH-A1FC-11A TCGA-E2-A158-11A
# RP11-368I23.2        0.5770281       0.00000000       0.33007291       0.93334894
# RP11-742D12.2        0.0000000       0.03483256       0.02550387       0.02258519
# RAB4B                3.0713037       3.90911352       3.50316793       2.76627551
# AC104183.2           0.0000000       0.00000000       0.00000000       0.00000000
expr=as.matrix(exp_filter)
#
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/得到模型后续分析/gse数据_hallmark通路富集评分_24_8_30.Rdata")
library(pheatmap)
pheatmap(gsva.kegg[1:50,], show_colnames = T, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# 设置分组
ann_col =Train_riskScore_cli[,c(1,14)] #创建分组列
rownames(ann_col)=ann_col$sample
ann_col=as.data.frame(ann_col[-1])
colnames(ann_col)="Sample"
row.names(ann_col) = colnames(expr) #这一行必须有，否则会报错：Error in check.length("fill") :  'gpar' element 'fill' must not be length 0

ann_color = list(Sample = c(High="#E95C59",Low="#4DBBD5E5")) #定义分组颜色

pheatmap(gsva.kegg[1:50,], show_colnames = T, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         annotation_col = ann_col, #表示是否对行、列进行注释，默认NA
         annotation = NA, annotation_colors = ann_color  #表示行注释及列注释的颜色，默认NA
         )
pheatmap(gsva.kegg[1:50,], 
         scale = "row", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         
         cluster_rows = T,cluster_cols = F, #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         
         cutree_rows = NA, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         
         treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         
         border_color = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         
         # cellwidth = 60, cellheight = 7.5,  #表示单个单元格的宽度\高度，默认为 “NA”
         
         display_numbers = F, #表示是否在单元格上显示原始数值或按照特殊条件进行区分标记
         
         fontsize_number = 6, #表示热图上显示数字的字体大小
         
         number_format = "%.2f", #表示热图单元格上显示的数据格式，“%.2f” 表示两位小数,“%.1e”表示科学计数法
         
         number_color = "grey30", #表示热图单元格上显示的数据字体颜色
         
         fontsize =10, fontsize_row = 6, fontsize_col = 10, #热图中字体大小、行、列名字体大小
         
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         
         main = "Gene标题", #表示热图的标题名字
         
         color = colorRampPalette(c("navy","white","firebrick3"))(100), #表示热图颜色,(100)表示100个等级
         
         angle_col = "45", #表示列标签的角度
         
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
         gaps_col = c(1,2,3,4,5,6),  #仅在未进行列聚类时使用，表示在列方向上热图的隔断位置
         
         # annotation_row = ann_row, 
         annotation_col = ann_col, #表示是否对行、列进行注释，默认NA
         
         annotation = NA, annotation_colors = ann_color,  #表示行注释及列注释的颜色，默认NA
         
         annotation_legend = TRUE, #表示是否显示注释的图例信息
         
         annotation_names_row = TRUE, annotation_names_col = TRUE) #表示是否显示行、列注释的名称

# 为啥分组都分开了，我觉得是顺序的问题
# 修改行注释的顺序
# 使用order函数对Sample列进行排序，使得"High"在前，"Low"在后
ann_col$group=NA # 必不可少，只有一列下面这个代码会报错
ann_col <- ann_col[order(ann_col$Sample, decreasing = TRUE), ]
ann_col=as.data.frame(ann_col[-2])
# 根据修改好得注释得顺序，修改gsva.kegg顺序
# 使用match函数找出mat1行名在mat2中的位置
match_indices <- match(rownames(ann_col), colnames(gsva.kegg)) #使用第一个矩阵的行名作为参照，找出这些行名在第二个矩阵中的位置
# 使用这个索引顺序来重新排列mat2
gsva.kegg <- gsva.kegg[,match_indices]
#
identical(rownames(ann_col), colnames(gsva.kegg))
#
pheatmap(gsva.kegg[1:50,], show_colnames = F, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         annotation_col = ann_col, #表示是否对行、列进行注释，默认NA
         annotation = NA, annotation_colors = ann_color  #表示行注释及列注释的颜色，默认NA
)
#
pheatmap(gsva.kegg[1:50,], show_colnames = F, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = F,
         color = colorRampPalette(c("#4DBBD5E5","white","#E95C59"))(100),
         annotation_col = ann_col, #表示是否对行、列进行注释，默认NA
         annotation = NA, annotation_colors = ann_color  #表示行注释及列注释的颜色，默认NA
)

#
save(ann_col,ann_color,gsva.kegg,file = "../出图_24_9_3/fig8d_hallmark_热图_画图需要.Rdata")


## 现在可以计算riskScore与通路得相关性了
rm(list=ls())
gc()
#
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/得到模型后续分析/gse数据_hallmark通路富集评分_24_8_30.Rdata")
#
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/不用封装的跑/选定模型计算好的riskscore_24_8_29.Rdata")
Train_riskScore_cli[1:4,1:4]
Train_riskScore_cli=Train_riskScore_cli[,c(1,13)]
rownames(Train_riskScore_cli)=Train_riskScore_cli$sample
Train_riskScore_cli=as.data.frame(Train_riskScore_cli[-1])
Train_riskScore_cli=t(Train_riskScore_cli)
Train_riskScore_cli=as.data.frame(Train_riskScore_cli)
#
identical(colnames(gsva.kegg),colnames(Train_riskScore_cli))
# [1] TRUE
exprSet=rbind(gsva.kegg,Train_riskScore_cli)
#
exprSet[49:51,1:4]
#                                     GSM1419942  GSM1419943 GSM1419944  GSM1419945
# HALLMARK_WNT_BETA_CATENIN_SIGNALING -0.2861151 -0.01692967  0.2692770 -0.03754189
# HALLMARK_XENOBIOTIC_METABOLISM       0.2257636 -0.01212682 -0.2393317  0.26050552
# riskScore                            1.9325508  1.47404423  3.1763459  1.04706084
exprSet=as.data.frame(exprSet)
class(exprSet[1,1])
# [1] "numeric"
#
exprSet1=exprSet
# 下面将行名去掉HALLMARK
rownames(exprSet)=str_split_fixed(rownames(exprSet),"_",n=2)[,2]
rownames(exprSet)[51]="riskScore"

# 2.写一个函数批量计算相关性
# 
# 这个函数只要输入一个基因，他就会批量计算这个基因跟其他编码基因的相关
# 
# 性，返回相关性系数和p值。

###对于有缺失值的基因，有效样本小于4会报错
batch_cor <- function(gene){
  y <- as.numeric(exprSet[gene,])
  rownames <- rownames(exprSet)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,type='spearman')
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
###这是修改的代码  加一个判断   样本量<10的就不要了吧
batch_cor <- function(gene){
  rownames <- rownames(exprSet)
  do.call(rbind,future_lapply(rownames, function(x){
    xy <- exprSet[c(gene,x),]
    xy <- t(xy) %>% na.omit() %>% as.data.frame()
    if (nrow(xy)>10){
      dd  <- cor.test(as.numeric(xy[,1]),as.numeric(xy[,2]),type='spearman')
      data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
    }
    
  }))
}

# 3.并行化运行函数
# 
# 以riskScore这个基因为例

library(future.apply)
# plan(multiprocess)
system.time(dd <- batch_cor('riskScore'))
# system.time(dd <- batch_cor('MFAP2'))
# 绘制热图
exprSet_t=t(exprSet)
corr <- cor(exprSet_t)
library(corrplot)
# 默认绘图样式
corrplot(corr)
col2 <- rev(COL2('RdBu', 100)) # 生成调色板后使用 rev() 函数来颠倒顺序
col2 = colorRampPalette(c("#4DBBD5E5","white","#E95C59"))(200)
corrplot(corr, method = c('pie'), 
         type = c('upper'), 
         col = col2, # 设置一个连续的
         outline = 'white', 
         # order = c('AOE'), 
         diag = TRUE,
         tl.cex = 0.5, #对角线文字大小
         tl.col = 'black', #对角线文字颜色
         tl.pos = 'td' # d仅在对角线显示文本标签
         # ,bg = "lightblue" # 设置图标背景颜色
         # ,mar = c(0,0,0,0)
)

save(exprSet,exprSet_t,corr,col2,file = "../出图_24_9_3/fig8e_hallmark相关性热图_表格.Rdata")

#### 计算hallmark上下调比较显著的通路与预后的关系。
rm(list=ls())
gc()
#
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/Mime1包/mydata/GSE58812_不分割log的清洗好的data.Rdata")
Train$OS.time=as.numeric(Train$OS.time)
Train=Train %>% filter(OS.time>200)
Train[1:4,1:4]
#
Train=Train[,c(1:3)]
head(Train)
#                    ID OS.time OS
# GSM1419942 GSM1419942    1520  1
# GSM1419943 GSM1419943    1281  1
# GSM1419944 GSM1419944    1066  1
# GSM1419945 GSM1419945    1050  1
# GSM1419946 GSM1419946     422  1
# GSM1419947 GSM1419947    2081  0
# 取hallmark的GSVA富集分析结果
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/得到模型后续分析/gse数据_hallmark通路富集评分_24_8_30.Rdata")
gsva=as.data.frame(t(gsva.kegg))
gsva[1:4,1:4]
#            HALLMARK_ADIPOGENESIS HALLMARK_ALLOGRAFT_REJECTION
# GSM1419942             0.3192342                   -0.3738940
# GSM1419943            -0.0299503                    0.0986041
# GSM1419944            -0.1203670                   -0.4587275
# GSM1419945             0.1270422                    0.2989670
identical(rownames(Train),rownames(gsva))
#
aimplot=cbind(Train,gsva)
#
aimplot$OS.time=as.numeric(aimplot$OS.time)
aimplot$OS=as.numeric(aimplot$OS)
#
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_ANGIOGENESIS>median(aimplot$HALLMARK_ANGIOGENESIS),'high','low')

aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION>median(aimplot$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_INFLAMMATORY_RESPONSE>median(aimplot$HALLMARK_INFLAMMATORY_RESPONSE),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_NOTCH_SIGNALING>median(aimplot$HALLMARK_NOTCH_SIGNALING),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_ALLOGRAFT_REJECTION>median(aimplot$HALLMARK_ALLOGRAFT_REJECTION),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_INTERFERON_ALPHA_RESPONSE>median(aimplot$HALLMARK_INTERFERON_ALPHA_RESPONSE),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_INTERFERON_GAMMA_RESPONSE>median(aimplot$HALLMARK_INTERFERON_GAMMA_RESPONSE),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_GLYCOLYSIS>median(aimplot$HALLMARK_GLYCOLYSIS),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_MYOGENESIS>median(aimplot$HALLMARK_MYOGENESIS),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_CHOLESTEROL_HOMEOSTASIS>median(aimplot$HALLMARK_CHOLESTEROL_HOMEOSTASIS),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_HYPOXIA>median(aimplot$HALLMARK_HYPOXIA),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_COAGULATION>median(aimplot$HALLMARK_COAGULATION),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_APICAL_JUNCTION>median(aimplot$HALLMARK_APICAL_JUNCTION),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_HEDGEHOG_SIGNALING>median(aimplot$HALLMARK_HEDGEHOG_SIGNALING),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_ESTROGEN_RESPONSE_EARLY>median(aimplot$HALLMARK_ESTROGEN_RESPONSE_EARLY),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_ESTROGEN_RESPONSE_LATE>median(aimplot$HALLMARK_ESTROGEN_RESPONSE_LATE),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_IL6_JAK_STAT3_SIGNALING>median(aimplot$HALLMARK_IL6_JAK_STAT3_SIGNALING),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_COMPLEMENT>median(aimplot$HALLMARK_COMPLEMENT),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_KRAS_SIGNALING_UP>median(aimplot$HALLMARK_KRAS_SIGNALING_UP),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_TNFA_SIGNALING_VIA_NFKB>median(aimplot$HALLMARK_TNFA_SIGNALING_VIA_NFKB),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_INFLAMMATORY_RESPONSE>median(aimplot$HALLMARK_INFLAMMATORY_RESPONSE),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_KRAS_SIGNALING_DN>median(aimplot$HALLMARK_KRAS_SIGNALING_DN),'high','low')
aimplot$Alveolar_quartile=ifelse(aimplot$HALLMARK_IL2_STAT5_SIGNALING>median(aimplot$HALLMARK_IL2_STAT5_SIGNALING),'high','low')

table(aimplot$Alveolar_quartile)
NK_OS <- survfit(Surv(OS.time,OS)~Alveolar_quartile,data = aimplot)
ggsurvplot(NK_OS,pval = T,risk.table = T,surv.median.line = 'hv',
           title='Overall survival',xlab='Days')

library(ggsci)
palette = "npg"  # 使用 Nature 期刊风格配色

ggsurvplot(NK_OS, data = aimplot,
           conf.int = TRUE,  pval = TRUE,
           surv.median.line = "hv",
           # risk.table = TRUE,
           # risk.table.height = 0.25,  # 调整风险表高度
           # risk.table.col = "strata", # 风险表按分层显示颜色
           palette = "npg",
           legend.labs=c("High","Low"), #标签
           legend.title="RiskScore",
           title="HALLMARK_IL2_STAT5_SIGNALING", #标题
           ylab="Overall survival",xlab = " Time (Days)", #更改横纵坐标
           censor.shape = 124,censor.size = 2, #删失点的形状和大小 break.x.by = 720#横坐标间隔 
           ggtheme = theme_minimal()+ # 使用干净的主题
                     theme(plot.title = element_text(hjust = 0.5),# 标题居中
                           panel.grid = element_blank(),  # 移除背景网格线
                           axis.line = element_line(color = "black", size = 0.5)  # 添加x轴和y轴线条
                           ) 
)