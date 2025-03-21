# 免疫浸润
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/得到模型后续分析/棒棒糖图展示风险评分与免疫浸润的关系_24_8_29.Rdata")
correlation[1:4,1:3] # 要下面这个表达矩阵即可
# 这个是前面计算相关性得到的，往前翻一下correlation就能找到了

data=as.data.frame(correlation)
colnames(data)=c("Cell","cor","pvalue")
data[1:4,1:3]

# 4 读取数据并根据数据确定图像数据
# 
# 4.1定义颜色

#定义圆圈颜色的函数
p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
} 

# 4.2定义棒棒糖圈大小

#定义设置圆圈大小的函数
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}
# 2.3根据定义的函数进行上色

#根据pvalue定义圆圈的颜色
points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color
points.cex = fcex(x=data$cor)
data$points.cex = points.cex

# 4.3.进行可视化

xlim = ceiling(max(abs(data$cor))*10)/10
pdf(file="../../../Breastcancer_EMBOJ/emboj_预后模型/101机器学习/得到模型后续分析/棒棒糖图_自己的数据_TNBC亚型数据.pdf", width=9, height=7)
pdf(file="../出图_24_9_3/fig7b_棒棒糖图_自己的数据_TNBC亚型数据.pdf", width=9, height=9)

# pdf(file="../Breastcancer_EMBOJ/new_analyse_24_4_17/ssgsea/棒棒糖图_自己的数据_TNBC亚型加了epi数据.pdf", width=9, height=7)
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2) 

# 直到这里我们就将基本的图形框架绘制完了，但是需要直观地展示仍需添加一些东西
# 添加图形的线段

#绘制图形的线段
segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
#绘制图形的圆圈
points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
#展示免疫细胞的名称
text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5) 

# 展示每个免疫细胞的p值

#展示pvalue
pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F) 

# 对图形添加图例

#绘制圆圈大小的图例
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off() 




library(CIBERSORT)
# 运行
f = "../ssgsea/ciber_GSE58812.Rdata" #设置一个默认的结果保存文件
if(!file.exists(f)){
  #devtools:: install_github ("Moonerss/CIBERSORT")
  library(CIBERSORT)
  lm22f = system.file("extdata", "LM22.txt", package = "CIBERSORT")
  TME.results = cibersort(lm22f, 
                          "../ssgsea/exp_cibersortUse.txt" ,  #存好的RNA-seq文件
                          perm = 1000, #迭代次数
                          QN = T) # # QN如果是芯片设置为T，如果是测序就设置为F，这边应该设置为F
  save(TME.results,file = f)
}
load(f)
TME.results[1:4,1:4]
#            B cells naive B cells memory Plasma cells T cells CD8
# GSM1419942    0.01064418              0   0.04222759  0.00000000
# GSM1419943    0.01363919              0   0.17402532  0.02386505
# GSM1419944    0.00000000              0   0.16042050  0.00000000
# GSM1419945    0.02107582              0   0.03571689  0.00000000

# 导出表格
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/ssgsea/ciber_GSE58812.Rdata")
write.csv(TME.results,file = "../出图_24_9_3/附表13_LM22_反卷积结果.csv")


# 添加风险分组
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/emboj_预后模型/101机器学习/不用封装的跑/选定模型计算好的riskscore_24_8_29.Rdata")
Train_riskScore_cli[1:4,1:4]
rownames(Train_riskScore_cli)=Train_riskScore_cli[,1]
#
identical(rownames(Train_riskScore_cli),rownames(TME.results))
# [1] TRUE
forinfiltrate=cbind(Train_riskScore_cli,TME.results)
colnames(forinfiltrate)

# colnames(forinfiltrate)[40]="EMT_like_CAF"
forinfiltrate[1:2,1:16]
#
df=forinfiltrate
# 将数据框转化为长格式，假设 riskScore2 是高低风险的分组标签
df_melted <- melt(df, id.vars = "riskScore2", 
                  measure.vars = c("B cells naive",               
                                   "B cells memory","Plasma cells","T cells CD8",                 
                                   "T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated",
                                   "T cells follicular helper","T cells regulatory (Tregs)","T cells gamma delta",         
                                   "NK cells resting","NK cells activated","Monocytes",                   
                                   "Macrophages M0","Macrophages M1","Macrophages M2",              
                                   "Dendritic cells resting","Dendritic cells activated","Mast cells resting",         
                                   "Mast cells activated","Eosinophils","Neutrophils"))

ggplot(df_melted, aes(x = variable, y = value, fill = riskScore2)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.6, position = position_dodge(width = 0.9)) + 
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) + 
  labs(x = "Cell Type", y = "Fraction", title = "") +
  scale_fill_manual(values = c("#E95C59","#4DBBD5E5"), name = "Risk Group") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  # theme_minimal()+
  theme_classic() +  # 设置为经典主题，没有背景
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 设置x轴标签垂直显示
    # axis.text.x = element_text(size = 10),  # 调整x轴标签字体大小
    axis.title.x = element_text(size = 12),  # 调整x轴标题字体大小
    axis.title.y = element_text(size = 12),  # 调整y轴标题字体大小
    panel.spacing = unit(1, "lines")  # 设置不同细胞类型之间的间隔
    ,panel.border = element_rect(colour = "gray", fill = NA, size = 1),  # 添加黑色边框
    axis.line.x = element_blank(),  # 去掉x轴线
    axis.line.y = element_blank(),  # 去掉y轴线
    axis.ticks.x = element_blank(),  # 去掉x轴刻度线
    axis.ticks.y = element_blank(),  # 去掉y轴刻度线
    axis.text.y = element_blank(),  # 去掉y轴刻度标签
    axis.title = element_blank()  # 去掉轴标题
  )



# TIDE小提琴图：
my_comparisons <- list( c("high", "low")) #添加比较分组
p1 <- ggviolin(tidy_res, x = 'riskGroup', y = 'TIDE', fill = 'riskGroup',
               palette = c("#E95C59","#4DBBD5E5"),
               width = 0.5,
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size=0.5, tip.length = 0.02, method = 't.test')+
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1),  # 设置x轴标签垂直显示
    # axis.text.x = element_text(size = 10),  # 调整x轴标签字体大小
    # axis.title.x = element_text(size = 12),  # 调整x轴标题字体大小
    # axis.title.y = element_text(size = 12),  # 调整y轴标题字体大小
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加黑色边框
    axis.line.x = element_blank(),  # 去掉x轴线
    axis.line.y = element_blank(),  # 去掉y轴线
    # axis.ticks.x = element_blank(),  # 去掉x轴刻度线
    # axis.ticks.y = element_blank(),  # 去掉y轴刻度线
    # axis.text.y = element_blank(),  # 去掉y轴刻度标签
    # axis.title = element_blank()  # 去掉轴标题
  ) 
p1




# 抗肿瘤免疫循环 雷达图
library(ggradar)
ggradar(leida, background.circle.transparency = 0, group.colours = c("#E95C59","#4DBBD5E5"),
        grid.min = -5,
        grid.mid = 1,
        grid.max = 5,  # 默认最大为1，如果画图数据表中有大于1的就报错了，得改
        values.radar = c(1,5,15)
)
