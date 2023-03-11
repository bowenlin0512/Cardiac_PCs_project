###SAN & public data analysis
library(harmony)
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(viridis) # scale_fill_viridis()
library(ggsignif)
library(viridis)
library(ggpubr)
library(Nebulosa)
###################
#*-------
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
fs=list.files(pattern = '.h5')
fs
sceList = lapply(fs, function(x){
  # x=fs[1]
  print(x)
  a=Read10X_h5( x )
  a[1:4,1:4] 
  library(stringr)
  (p=str_split(x,'_',simplify = T)[,2])
  sce <- CreateSeuratObject( a ,project = p )
  sce
})
#*-------
###############################################merge####################################################################
SAN<-readRDS("scRNA_harmony_before_definaiton.rds")
DefaultAssay(immune.combined) <- "integrated"
DimPlot(SAN, label = T)
a$orig.ident<-"E14.5"
a$sample<-"E14.5"
a$tech<-"E14.5"
scRNA_harmony <- merge(SAN,a)
#DefaultAssay(scRNA_harmony) <- "integrated"
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "tech")})
#system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
ElbowPlot(scRNA_harmony, ndims = 50)
pc.num=1:40
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:40)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:40) %>% FindClusters()
#group_by_cluster
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T) 
#group_by_sample
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident') 
#combinate
plotc <- plot1+plot2
ggsave("scRNA_harmony.png", plot = plotc, width = 10, height = 5)
#saveRDS(scRNA_harmony, 'scRNA_harmony_merge_E13.5.rds')
scRNA_harmony<-readRDS("scRNA_harmony_merge_E13.5.rds")
DimPlot(scRNA_harmony, reduction = "umap", label=T) 
DimPlot(scRNA_harmony, reduction = "umap", split.by = "sample",label=T) 
###findmarker
markers.to.plot <- c("T","Tbx6","Kdr","Mixl1","Eomes","Gsc",
                     "Mesp1","Tbx3","Tbx5","Tbx18","Isl1",
                     "Osr1","Col1a1","Col3a1","Myl2","Myl1","Myl7","Tnnt2",
                     "Myh6","Vsnl1","Shox2","Hcn4","Cacna2d2","Afp","Klf2")
FeaturePlot(scRNA_harmony, features = markers.to.plot,label = TRUE)
FeaturePlot(scRNA_harmony, features = markers.to.plot,split.by = "sample",label = F)
p<-FeaturePlot(scRNA_harmony, features = c("Hcn4","Shox2"), split.by = "sample",label = TRUE)
p+scale_color_viridis()
markers.to.plot <- c("Hcn4","Shox2")
##cividis,magma,inferno,viridis
DotPlot(scRNA_harmony, features = rev(markers.to.plot), dot.scale = 8) + RotatedAxis()+scale_color_viridis(option="cividis")
DotPlot(scRNA_harmony, features = rev(markers.to.plot), dot.scale = 8) + scale_color_viridis()
DotPlot(scRNA_harmony, features = rev(markers.to.plot)) + scale_color_viridis()
DotPlot(SAN, features = rev(markers.to.plot)) + scale_color_viridis()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DimPlot(SAN, reduction = "umap", label=T) 
VlnPlot(scRNA_harmony, features = markers.to.plot,log = F,slot = "data",pt.size = 0)
VlnPlot(scRNA_harmony, features = markers.to.plot,log = F,slot = "data",pt.size = 0.1)
seurat.object_copy <- scRNA_harmony
my_levels <- c("d5","d8","d13","E14.5")
seurat.object_copy$sample <- factor(x = seurat.object_copy$sample, levels = my_levels)
FeaturePlot(seurat.object_copy, features = markers.to.plot,split.by = "sample",label = TRUE)
FeaturePlot(seurat.object_copy, features = markers.to.plot,split.by = "sample",label = F)
scRNA_harmony.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scRNA_harmony.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(scRNA_harmony.markers,"scRNA_harmony.markers.csv")
top20 <- scRNA_harmony.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(scRNA_harmony, features = top20$gene) + NoLegend()
######################################################################################################
#####Figrues
Cells.sub <- subset(seurat.object_copy@meta.data, sample=="d5"|sample=="d8"|sample=="d13")
SAN_sub <- subset(seurat.object_copy, cells=row.names(Cells.sub))
DimPlot(SAN, reduction = "umap", label = T, split.by = "sample",pt.size = 0.5)
markers.to.plot <- c("T","Tbx6","Kdr","Mixl1","Eomes","Gsc",
                     "Mesp1","Tbx3","Tbx5","Tbx18","Isl1",
                     "Osr1","Col1a1","Col3a1","Myl2","Myl1","Myl7","Tnnt2",
                     "Myh6","Vsnl1","Shox2","Hcn4","Cacna2d2","Afp","Klf2",
                     "Gata4","Gata5","Gata6","Pdgfra","Rgs6","Hand2",
                     "Myl4","Nppa","Myh9","Cnn1","Dusp27","Fhl2","Fhl1",
                     "Sox5", "Pcdh7", "Prelid2", "Grxcr2","Irx4","Mki67","Foxa2",
                     "Cdh5","Mfap2")
DotPlot(SAN, features = rev(markers.to.plot)) + scale_color_viridis()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
FeaturePlot(SAN_sub, features = markers.to.plot,split.by = "sample",label = F)
##Change the colors
library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1"), 
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00",
                    "#7CFC01","#7CFC02")  
DimPlot(SAN_sub, label = T, pt.size = 1,cols = cell_type_cols)+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
library(ggsci)
cell_type_cols <- c(pal_npg("nrc", alpha = 0.6)(9), 
                    pal_jama("default", alpha = 0.6)(7),
                    pal_jco("default", alpha = 0.6)(9),
                    pal_d3("category20")(20))  
cell_type_cols <- c(pal_d3("category20")(16))  
cell_type_cols <- c(pal_jco("default", alpha = 0.6)(7),
                    pal_npg("nrc", alpha = 0.6)(9))  
DimPlot(SAN_sub, reduction = "umap", label = T, pt.size = 0.5,cols = cell_type_cols)
DimPlot(SAN_sub, reduction = "umap", label = T, pt.size = 0.5)

##Check the genes
DimPlot(SAN, reduction = "umap", label = T, split.by = "sample",pt.size = 0.5)
markers.to.plot <- c("T","Tbx6","Hoxb1","Snai1","Kdr","Mixl1","Lhx1","Eomes","Gsc",
                     "Mesp1","Tbx3","Tbx5","Tbx18","Isl1","Col1a1","Col3a1","Myl2","Myl1","Myl7","Tnnt2",
                     "Myh6","Vsnl1","Shox2","Hcn4","Cacna2d2","Afp","Klf2",
                     "Gata4","Gata5","Gata6","Pdgfra","Rgs6","Hand2",
                     "Myl4","Myh9","Cnn1","Dusp27","Fhl2","Fhl1",
                     "Pcdh7", "Prelid2", "Irx4","Mki67","Epcam","Foxa2",
                     "Cldn7","Cdh5","Mfap2")

# list_genes=split(topN$gene, topN$cluster)
list_genes=list(T=c("PTPRC","CD3D","CD3E", "CD3G", "CD4", "CD8A","CD8B",
                    "CCR7", "LEF1", "GZMH", "GZMK"),
                B=c("CD79A", "CD79B", "MS4A1", "TCL1A", "CD22", "CD19"),
                mono=c("CD14",  "S100A8", "S100A9", "FCN1",'LYZ',
                       "FCGR3A", "MS4A7" ),
                NK=c("NKG7","GZMB", "GZMA", "CST7"),
                DC=c("FCER1A", "CST3", "CLEC10A"),
                platelet=c("PPBP","GP9","ITGA2B","PF4" ) )


DotPlot(SAN, group.by = 'cell_type',features = rev(markers.to.plot)) + scale_color_viridis()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##dotplot
#devtools::install_github("Simon-Leonard/FlexDotPlot")
library(FlexDotPlot)
dp=DotPlot(SAN, features = markers.to.plot,group.by = 'cell_type') + RotatedAxis()
# Dot plot with shape type (and not size) controlled by "Percent Expressed"
dot_plot(dp$data[,c(3,4,1,2,5)], size_var = "pct.exp", col_var = "avg.exp.scaled",
         size_legend = "Percent Expressed",col_legend = "Average expression",
         x.lab.pos = "bottom",y.lab.pos = "right",display_max_sizes = F,
         shape.scale = 8, hclust_method = "ward.D2",
         dend_x_var = c("pct.exp","avg.exp.scaled"),
         dend_y_var = c("pct.exp","avg.exp.scaled"),
         text.size = 0.5,
         text.vjust = 0.5,
         size.breaks.number = 4,
         color.breaks.number = 3,
         #cols.use = c("#330066","#336699","#66CC66","#FFCC33"),
         cols.use = c("#4DBBD599","#FFDC9199","#CC3333"),
         x.lab.size.factor = 1.5,
         y.lab.size.factor = 1,
         transpose=F)
DotPlot(SAN, features = rev(markers.to.plot)) + scale_color_viridis()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeaturePlot(SAN_sub, features = markers.to.plot,split.by = "sample",label = F)
cluster2celltype <- c("0"="Nascent mesoderm", 
                      "1"="Nascent mesoderm", 
                      "2"="Ventricular-like myocyte", 
                      "3"= "Transitional cell", 
                      "4"= "Atrium-like myocyte", 
                      "5"= "Ventricular-like myocyte",
                      "6"= "Transitional cell", 
                      "7"= "Transitional cell", 
                      "8"= "SAN",
                      "9"="Mesendoderm/epiblast",
                      "10"="Pacemaker progenitor cell", 
                      "11"="Primitive streak-derived tissues", 
                      "12"= "Nascent mesoderm", 
                      "13"= "Cardiac fibroblast", 
                      "14"= "Immature cardiomyocyte",
                      "15"= "Fibroblast", 
                      "16"= "Cardiac progenitor cell", 
                      "17"= "Mesendoderm/epiblast",
                      "18"="Cardiac fibroblast",
                      "19"="Paraxial mesoderm", 
                      "20"="SHF", 
                      "21"= "Endoderm", 
                      "22"= "Endothelial")
SAN[['cell_type']] = unname(cluster2celltype[SAN@meta.data$seurat_clusters])
DimPlot(SAN, reduction = 'umap', group.by = 'cell_type',
        label = TRUE, pt.size = 0.5,cols = cell_type_cols,repel = T)
DimPlot(SAN, reduction = 'umap', group.by = 'cell_type',
        label = TRUE, pt.size = 1.5,cols = cell_type_cols,repel = T)
DimPlot(SAN, reduction = 'umap', group.by = 'cell_type',
        label = TRUE, pt.size = 1,cols = cell_type_cols,repel = T)
DimPlot(SAN, reduction = 'umap', group.by = 'cell_type',split.by = 'sample',
        label = T, pt.size = 1,cols = cell_type_cols,repel = T)

##
seurat.object_copy <- SAN
my_levels <- c("d5","d8","d13")
seurat.object_copy$sample <- factor(x = seurat.object_copy$sample, levels = my_levels)
DimPlot(seurat.object_copy, reduction = 'umap', split.by = 'sample',
        label = F, pt.size = 0.5,repel = T)
DimPlot(seurat.object_copy, reduction = 'umap', group.by = 'sample',label = F, pt.size = 0.5,cols = cell_type_cols,repel = T)
DimPlot(seurat.object_copy, reduction = 'umap', group.by = 'cell_type',split.by = 'sample',cols = cell_type_cols,
        label = F, pt.size = 0.5,repel = T)
DimPlot(seurat.object_copy, reduction = 'umap', group.by = 'cell_type',split.by = 'tech',cols = cell_type_cols,
        label = F, pt.size = 0.5,repel = T)
FeaturePlot(seurat.object_copy, features = markers.to.plot,split.by = "sample",label = TRUE)
cols<-brewer.pal(5, 'Blues')
FeaturePlot(seurat.object_copy, features = markers.to.plot,label = TRUE,cols = cols)
FeaturePlot(seurat.object_copy, features = markers.to.plot,label = TRUE,pt.size = 0.5,cols = c("grey","red"))
##
ggplot(seurat.object_copy@meta.data,aes(x=sample,fill=cell_type))+
  geom_bar(position = "fill",width = 0.7)+
  xlab(NULL)+ylab("cell.percentage")+
  scale_fill_manual(values =cell_type_cols)+theme_bw()

##
markers.to.plot <- c("Vsnl1","Bmp4","Socs2","Dsp")
markers.to.plot <- c("Zfp503")
FeaturePlot(seurat.object_copy, features = markers.to.plot,  split.by = "sample", label = TRUE,repel = T)
FeaturePlot(seurat.object_copy, features = markers.to.plot,label = F,repel = T)
FeaturePlot(seurat.object_copy, features = markers.to.plot, split.by = "sample", cols = viridis(20),label = TRUE,repel = T)
FeaturePlot(seurat.object_copy, features = markers.to.plot,label = TRUE,repel = T)
FeaturePlot(object = seurat.object_copy, features = markers.to.plot,
            cols = c("lightgrey" ,"blue"),slot = "data",
            label.size = 6,pt.size = 1.2) + 
  theme(axis.line = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.title = element_blank()) + 
  guides(color=F)

##
plot_density(immune, "CD4")
plot_density(seurat.object_copy, c("Hcn4", "Shox2"),joint = TRUE,pal = "viridis")
plot_density(seurat.object_copy, c("Zfp503"),joint = TRUE,pal = "viridis")
plot_density(seurat.object_copy, c("Cacna2d2", "Igfbp5","Lbh", "Dact1"),joint = F,pal = "viridis")
plot_density(seurat.object_copy, c("Hcn4", "Shox2"),joint = T,pal = "viridis")
##Violin plots
VlnPlot(seurat.object_copy, features = "Shox2",group.by = 'sample',pt.size = 0.5)+NoLegend()
VlnPlot(seurat.object_copy, features = "Shox2")
VlnPlot(seurat.object_copy, features = markers.to.plot,group.by = 'cell_type',pt.size = 0,cols = cell_type_cols,ncol =2)+NoLegend()
VlnPlot(seurat.object_copy, features = "Zfp503",group.by = 'cell_type',pt.size = 0,cols = cell_type_cols,ncol =2)+NoLegend()
VlnPlot(seurat.object_copy, features = c("Hcn4", "Shox2","Vsnl1"),group.by = 'sample',pt.size=0.2)+
  NoLegend()+scale_fill_npg()
p<-VlnPlot(seurat.object_copy, features = "Hcn4",group.by = 'sample',pt.size=0.2)
DotPlot(seurat.object_copy, group.by = 'cell_type',features = rev(markers.to.plot)) + scale_color_viridis()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####
Vlnpubr <- function(seo, gene_signature, file_name, test_sign,group_my,label  = "p.format"){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(seo , features = signature,
            pt.size = 0.1, 
            group.by =  group_my, 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) +  stat_compare_means(comparisons = test_sign, label = label ) + stat_compare_means()+NoLegend()
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  #file_name <- paste0(file_name, "_r.png")
  #ggsave(file_name, width = 14, height = 8)
}
##Compare the expression
gene_sig <- c("Hcn4","Shox2","Ryr2")
comparisons <- list(c("d5", "d8"),c("d5", "d13"),c("d8", "d13"))
Vlnpubr(seurat.object_copy,gene_signature = gene_sig, test_sign = comparisons,group_my = 'sample')
#############
#saveRDS(seurat.object_copy, file = "scRNA_harmony_order_right.rds")
seurat.object_copy<-readRDS("scRNA_harmony_order_right.rds")