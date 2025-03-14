# scRNA
sce <- readRDS("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/new_analyse_24_4_17/sce_primary_er_her2_tnbc_normal_lastFilter_24_4_17.rds")
sce <- NormalizeData(sce, normalization.method =  "LogNormalize", 
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", nfeatures = 3000) 
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
res.used <- 1.0
ElbowPlot(sce,ndims = 50)
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(object = sce, verbose = T, resolution = res.used)
# sce <- FindClusters(object = sce, verbose = T, resolution = 0.5)
set.seed(123)
sce <- RunUMAP(object = sce, dims = 1:20, do.fast = TRUE)
sce <- RunTSNE(object = sce, dims = 1:20, do.fast = TRUE)
# DimPlot(sce,reduction = "umap",label=T,raster=FALSE)
# DimPlot(sce,reduction = "umap",label=T, group.by = "orig.ident")
table(sce@meta.data$seurat_clusters) 

p1=DimPlot(sce,group.by = "seurat_clusters",label = T,raster=FALSE)
p2=DimPlot(sce,group.by = "orig.ident",label = T,raster=FALSE)+NoLegend()
p3=DimPlot(sce,group.by = "Tissue",label = T,raster=FALSE)
p1+p2+p3

p1=DimPlot(sce,group.by = "seurat_clusters",label = T,raster=FALSE,reduction = "tsne")+NoLegend()
p2=DimPlot(sce,group.by = "orig.ident",label = T,raster=FALSE,reduction = "tsne")+NoLegend()
p3=DimPlot(sce,group.by = "Tissue",label = T,raster=FALSE,reduction = "tsne")
p1+p2+p3


# 注释看看
all_nonimmune=list(
  immune=c("PTPRC"), 
  tcell=c("CD3D","CD3E"),
  nk=c("NCAM1","GNLY","NKG7","KLRD1"),
  Stromal=c("DCN","COL1A2","COL1A1","LUM","PDGFRA","PDGFRB"
            ),
  Plasma_bcell=c("MZB1","IGHG4"),
  Myepithelials=c("KRT14","ACTA2"),
  myeloid=c("CD68","CD163","LYZ","SPP1","CST3",
            "LST1","LILR82", 
            "C1QC","C1QA", 
            "TREM2"), 
  epi=c("CD24","EPCAM","KRT19","KRT7","KRT8","KRT18"
        # ,"KRT5"
  ), 
  edo=c("PECAM1","VWF","ENG","CDH5","PLVAP"), # 内皮
  CD20_bcell=c("MS4A1","CD79A")
)
th=theme(axis.text.x = element_text(angle = 90, 
                                    vjust = 0.5, hjust=0.5))  

p <- DotPlot(sce, features = all_nonimmune,
             assay='RNA')+  th

p
# 注释一下
celltype=data.frame(ClusterID=0:48,
                    celltype="Epithelials") 
celltype[celltype$ClusterID %in% c(3,32),2]='Edothelials'  #
celltype[celltype$ClusterID %in% c(0,8,39,41),2]='T_NK'  # 
celltype[celltype$ClusterID %in% c(23),2]='Naive_Bcells' # 

celltype[celltype$ClusterID %in% c(22,48),2]='Plasma_Bcells' #

celltype[celltype$ClusterID %in% c(14,16,38),2]='Myeloids' 
celltype[celltype$ClusterID %in% c(40),2]='Mast' 
celltype[celltype$ClusterID %in% c(6,9,13,30,46),2]='Fibroblasts'
celltype[celltype$ClusterID %in% c(34),2]='Myepithelials' # 即epi+fib
celltype[celltype$ClusterID %in% c(19),2]='Undefined'
celltype[celltype$ClusterID %in% c(1,2,4,5,7,10,11,12,15,17,18,20,21,24,25,26,27,28,29,
                                   31,33,35,36,37,42,43,44,45,47),2]='Epithelials' 

head(celltype)
celltype 
table(celltype$celltype)
sce@meta.data$firstannotation = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'firstannotation'] <- celltype$celltype[i]}
table(sce@meta.data$firstannotation)

DimPlot(sce, reduction = "tsne",cols = mycolor2,group.by = "firstannotation", label = TRUE,pt.size = 0.5,label.box = T,raster = FALSE)




# plot1cell
library(plot1cell)
#
sce <- readRDS("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/new_analyse_24_4_17/sce_primary_er_her2_tnbc_normal_第一次注释_24_4_18.rds")
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/new_analyse_24_4_17/phe_大群与T_NK_myeloids_Bcell整合_24_6_17.Rdata")
sce@meta.data=merge_meta
# 
head(sce@meta.data)
#
DimPlot(sce, reduction = "tsne",group.by = "firstannotation_modify", label = TRUE,pt.size = 0.5,label.box = T,raster = FALSE)
DimPlot(sce, reduction = "umap",group.by = "firstannotation_modify", label = TRUE,pt.size = 0.5,label.box = T,raster = FALSE)
DimPlot(sce, reduction = "umap",group.by = "Tissue", label = TRUE,pt.size = 0.5,label.box = T,raster = FALSE)

#
head(sce,2)
Idents(sce)=sce$firstannotation_modify


###Check and see the meta data info on your Seurat object
colnames(sce@meta.data)  

# turn tsne/umap
source("../出图_24_9_3/plot_circlize.R")

###Prepare data for ploting
circ_data <- prepare_circlize_data(sce, scale = 0.8 )
colnames(circ_data)
#
set.seed(1234)
#
#设置需要的颜色
mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#7F2268','#91D0C1','#bebcdf')
mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#7F2268','#91D0C1','#FFCC4F')

cluster_colors<-mycolors
group_colors<-rand_color(length(names(table(sce$Tissue)))) #分组
# "#145779FF" "#8A04D0FF" "#467F08FF" "#B6D028FF"
# "#A008CAFF" "#CD590FFF" "#DF83D8FF" "#DFC5FBFF"
rep_colors<-rand_color(length(names(table(sce$orig.ident))))

###plot and save figures

plot_circlize(circ_data,do.label = T, pt.size = 0.02,
              col.use = cluster_colors ,bg.color = '#f9f9f9e5', #'#f8f2e4'
              kde2d.n = 200, repel = T, label.cex = 1.2)

add_track(circ_data, group = "Tissue", colors = group_colors, track_num = 2) 

add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) 

# 修改plot_circlize.R脚本可以达到修改字体大小，以及对tsne进行画图，默认是umap

# 对tsne结果进行可视化



#  添加亚群可视化结果
sce_fib <- readRDS("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/new_analyse_24_4_17/sce_primary_er_her2_tnbc_normal_fib_24_4_19.rds")
DimPlot(sce_fib, reduction = "tsne",group.by = "fibannotation", label = TRUE,pt.size = 0.5,label.box = T,raster = FALSE)
DimPlot(sce_fib, reduction = "umap",group.by = "fibannotation", label = TRUE,pt.size = 0.5,label.box = T,raster = FALSE)

sce_bcell <- readRDS("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/new_analyse_24_4_17/sce_bcell_all_24_6_16.rds")
#
DimPlot(sce_bcell, reduction = "tsne",group.by = "bcellannotation", label = TRUE,pt.size = 0.5,label.box = T,raster = FALSE)
DimPlot(sce_bcell, reduction = "umap",group.by = "bcellannotation", label = TRUE,pt.size = 0.5,label.box = T,raster = FALSE)

sce_NK <- readRDS("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/new_analyse_24_4_17/sce_从T里面分出NK_24_6_16.rds")
# 
DimPlot(sce_NK, reduction = "tsne",group.by = "t_nkcellannotation", label = TRUE,pt.size = 0.5,raster = FALSE)
# 这边需要把T细胞注释的结果映射到T_NK注释中
load("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/new_analyse_24_4_17/phe_T_NK整合_24_6_17.Rdata")
sce_NK@meta.data=merge_meta
DimPlot(sce_NK, reduction = "tsne",group.by = "t_nk_integrated_6_17", label = TRUE,pt.size = 0.5,raster = FALSE)
DimPlot(sce_NK, reduction = "umap",group.by = "t_nk_integrated_6_17", label = TRUE,pt.size = 0.5,raster = FALSE)
#
sce_myeloids <- readRDS("~/my_prepare/lfr_finally/Breastcancer_EMBOJ/new_analyse_24_4_17/sce_myeloids_24_9_5.rds")
DimPlot(sce_myeloids, reduction = "tsne",group.by = "myeloidannotation", label = TRUE,pt.size = 0.5,raster = FALSE)
DimPlot(sce_myeloids, reduction = "umap",group.by = "myeloidannotation", label = TRUE,pt.size = 0.5,raster = FALSE)

# 组成list
sub.celltype_list=list(sce_fib,sce_bcell,sce_NK,sce_myeloids)
names(sub.celltype_list)=c("Fib","Bcell","T_NK","Myeloids")
save(sub.celltype_list,file = '../出图_24_9_3/fig1_绘图需要.Rdata')
# 添加四周亚型umap图
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
###Fibroblast subtypes
Fibroblast <- sub.celltype_list$Fib
Idents(Fibroblast) <- "fibannotation"
subcolors <- my36colors[1:nlevels(Fibroblast)]
#subcolors <- c('#bff542','#83f78f','#EBA1A2','#D70016','#eab3fc','#83b1f7','#D70016','#eab3fc','#83b1f7')
Fibroblast_meta<-get_metadata(Fibroblast, color = subcolors)

Fibroblast_meta %>%
  dplyr::group_by(fibannotation) %>%
  summarise(x = median(x),y = median(y)) -> centers_Fib

points(Fibroblast_meta$x*0.32-1.3,Fibroblast_meta$y*0.32-0.74, pch = 19, col = alpha(Fibroblast_meta$Colors,0.5), cex = 0.1); # x减的越小，越靠右边,y减的越小，越靠上
text(centers_Fib$x*0.32-1.3,centers_Fib$y*0.32-0.74, labels=centers_Fib$fibannotation, 
     cex = 0.8 # 字体大小
     , col = 'black')

# 注意这里的subcolors 可以自定义，也可以每次都使用my36colors 中的颜色，但是一定要注意以下2点


#T subtypes
T.sub <- sub.celltype_list$T_NK
Idents(T.sub) <- "t_nk_integrated_6_17"
subcolors <- my36colors[1:nlevels(T.sub)]
T_meta<-get_metadata(T.sub, color = subcolors)
T_meta %>%
  dplyr::group_by(t_nk_integrated_6_17) %>%
  summarise(x = median(x = x),y = median(x = y)) -> centers_T
points(T_meta$x*0.32+1.2,T_meta$y*0.32+0.73, pch = 19, col = alpha(T_meta$Colors,0.5), cex = 0.1);
text(centers_T$x*0.32+1.2,centers_T$y*0.32+0.73, labels=centers_T$t_nk_integrated_6_17, cex = 0.6, col = 'black')

#Myeloid subtypes
Myeloid.sub <- sub.celltype_list$Myeloids
Idents(Myeloid.sub)="myeloidannotation"
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
subcolors <- my36colors[1:nlevels(Myeloid.sub)]
# subcolors <- c('#D6E7A3','#E5D2DD', '#53A85F', '#F1BB72','#57C3F3', '#F3B1A0', '#476D87')
Myeloid_meta<-get_metadata(Myeloid.sub, color = subcolors)
Myeloid_meta %>%
  dplyr::group_by(myeloidannotation) %>%
  summarise(x = median(x = x),y = median(x = y)) -> centers_Mye
points(Myeloid_meta$x*0.32-1.2,Myeloid_meta$y*0.32+0.73, pch = 19, col = alpha(Myeloid_meta$Colors,0.5), cex = 0.1);
text(centers_Mye$x*0.32-1.2,centers_Mye$y*0.32+0.73, labels=centers_Mye$myeloidannotation, cex = 0.6, col = 'black')

##B subtype
Epi.sub <- sub.celltype_list$Bcell
Idents(Epi.sub)="bcellannotation"
subcolors <- my36colors[1:nlevels(Epi.sub)]
Epi_meta<-get_metadata(Epi.sub, color = subcolors)
Epi_meta %>%
  dplyr::group_by(bcellannotation) %>%
  summarise(x = median(x = x),y = median(x = y)) -> centers_Epi

points(Epi_meta$x*0.3+1.2,Epi_meta$y*0.3-0.73, pch = 19, col = alpha(Epi_meta$Colors,0.5), cex = 0.1);
text(centers_Epi$x*0.3+1.2,centers_Epi$y*0.3-0.73, labels=centers_Epi$bcellannotation, cex = 0.6, col = 'black')

#
# 3 ，添加四周umap的title 和 track的legend 
# 
# （1）添加，优化四周umap的title ，注意位置和大小

title_text <- function(x0, y0, x1, y1, text, rectArgs = NULL, textArgs = NULL) {
  center <- c(mean(c(x0, x1)), mean(c(y0, y1)))
  do.call('rect', c(list(xleft = x0, ybottom = y0, xright = x1, ytop = y1), rectArgs))
  do.call('text', c(list(x = center[1], y = center[2], labels = text), textArgs))
}

title_text(x0 = -1.35, x1 = -1.05, y0 = -1.06, y1=-1, text = 'Stromal',
           rectArgs = list(border='#F9F2E4',lwd=0.5),
           textArgs = list(col='black',cex = 1))

title_text(x0 = 1.05, x1 = 1.35, y0 = -1.06, y1=-1, text = 'Bcells',
           rectArgs = list(border='#F9F2E4',lwd=0.5),
           textArgs = list(col='black',cex = 1))

title_text(x0 = -1.35, x1 = -1.05, y0 = 1.06, y1=1, text = 'Myeloids',
           rectArgs = list(border='#F9F2E4',lwd=0.5),
           textArgs = list(col='black',cex = 1))

title_text(x0 = 1.05, x1 = 1.35, y0 = 1.06, y1=1, text = 'T_NK cells',
           rectArgs = list(border='#F9F2E4',lwd=0.5),
           textArgs = list(col='black',cex = 1))
# （2）添加track的legend

#plot group#
col_use<-c('#00288A','#DD001F','#84D000','#00CB47','#947F00','#006234')
col_use=c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
          '#7F2268','#91D0C1','#FFCC4F')
# mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
#              '#7F2268','#91D0C1','#FFCC4F')
cc<-get_metadata(sce, color = col_use)
cc %>%
  dplyr::group_by(Tissue) %>%
  summarise(x = median(x = x),y = median(x = y)) -> centers
# group_colors<-rand_color(length(names(table(sce$Tissue)))) #分组 ,颜色前面搞过了
col_group<-c("#145779FF","#8A04D0FF","#467F08FF","#B6D028FF")
lgd_points = Legend(labels = names(table(cc$Tissue)), type = "points", 
                    title_position = "topleft", 
                    title = "Group",
                    title_gp = gpar(col='black',fontsize = 7, fontface='bold'),
                    legend_gp = gpar(col = col_group),
                    labels_gp = gpar(col='black',fontsize = 5),
                    grid_height = unit(2, "mm"),
                    grid_width = unit(2, "mm"),
                    background = col_group)
draw(lgd_points, x = unit(15, "mm"), y = unit(50, "mm"),
     just = c("right", "bottom"))