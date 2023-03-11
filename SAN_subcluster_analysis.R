##scRNA-seq further analysis
##提取d5 d8的sample
seurat.object_copy@meta.data$celltype<-scRNA_harmony@active.ident
Cells.sub <- subset(scRNA_harmony@meta.data, sample=="d5"|sample=="d8")
seurat.object_copy <- subset(scRNA_harmony, cells=row.names(Cells.sub))

#提取基因表达
select_cells <- WhichCells(seurat.object_copy, slot = 'counts', expression = Shox2 > 0 | Hcn4 > 0) 
select_obj <- subset(seurat.object_copy, cells = select_cells)

doble_pos=colnames(subset(x = seurat.object_copy, subset = Shox2 > 0 & Hcn4 > 0, slot = 'counts'))
select_obj@meta.data$state=ifelse(colnames(select_obj) %in% doble_pos,'doble_pos','sigle_pos')
  
table(select_obj@meta.data$state)
markers <- FindMarkers(select_obj, ident.1 = "high", 
                       group.by = 'highORlow', 
                       subset.ident = "0")

markers <- FindMarkers(select_obj ,  ident.1="doble_pos", ident.2="sigle_pos",
                       group.by = 'state',logfc.threshold = 0.25)
write.csv(markers,"doble_pos_vs_sigle_pos.csv")
head(x = markers)


##提取Tbx18 Hcn4+细胞
#提取基因表达
select_cells <- WhichCells(seurat.object_copy, slot = 'counts', expression = Shox2 > 0 & Hcn4 > 0 ) 
select_obj <- subset(seurat.object_copy, cells = select_cells)

doble_pos=colnames(subset(x = seurat.object_copy, subset = Tbx18 > 0 & Hcn4 > 0, slot = 'counts'))
select_obj@meta.data$state=ifelse(colnames(select_obj) %in% doble_pos,'doble_pos','sigle_pos')

table(select_obj@meta.data$state)
markers <- FindMarkers(select_obj, ident.1 = "high", 
                       group.by = 'highORlow', 
                       subset.ident = "0")

markers <- FindMarkers(select_obj ,  ident.1="doble_pos", ident.2="sigle_pos",
                       group.by = 'state',logfc.threshold = 0.25)
write.csv(markers,"Tbx18_doble_pos_vs_sigle_pos.csv")
head(x = markers)

##compare 2
##hcn+
select_cells <- WhichCells(seurat.object_copy, slot = 'counts', expression = Hcn4 > 0 & Shox2 > 0) 
select_obj <- subset(seurat.object_copy, cells = select_cells)

doble_pos=colnames(subset(x = select_obj, subset = Tbx18 > 0 , slot = 'counts'))
select_obj@meta.data$state=ifelse(colnames(select_obj) %in% doble_pos,'doble_pos','sigle_pos')

table(select_obj@meta.data$state)
markers <- FindMarkers(select_obj, ident.1 = "high", 
                       group.by = 'highORlow', 
                       subset.ident = "0")

markers <- FindMarkers(select_obj ,  ident.1="doble_pos", ident.2="sigle_pos",
                       group.by = 'state',logfc.threshold = 0.25)
write.csv(markers,"Hcn4Shox2+_Tbx18-pos_vs_Tbx18-neg.csv")
head(x = markers)




##提取d13 Shox2及Hcn4 细胞
double_pos=colnames(subset(x = seurat.object_copy, subset = Tbx18 > 0  & Hcn4 > 0 , slot = 'counts'))
Tbx18_pos=colnames(subset(x = seurat.object_copy, subset = Tbx18 > 0 & Hcn4 <= 0 , slot = 'counts'))
Hcn4_pos=colnames(subset(x = seurat.object_copy, subset = Tbx18 <= 0 & Hcn4 > 0 , slot = 'counts'))

highORlow=ifelse(colnames(seurat.object_copy) %in% double_pos,'double_pos',
                 ifelse(colnames(seurat.object_copy) %in% Tbx18_pos,'Tbx18_pos',
                        ifelse(colnames(seurat.object_copy) %in% Hcn4_pos,'Hcn4_pos','double_neg')))
table(highORlow)
seurat.object_copy@meta.data$highORlow=highORlow
##去除doube_neg的细胞
Cells.sub <- subset(seurat.object_copy@meta.data, highORlow=="double_pos"|highORlow=="Tbx18_pos"|highORlow=='Hcn4_pos')
seurat.object_copy2 <- subset(seurat.object_copy, cells=row.names(Cells.sub))

##直接输入相关基因
cg_markers<-c("Hcn4","Tbx18","Shox2","Tbx3","Vsnl1","Isl1","Tbx5","Cd166", ##起搏
              "Myl3","Actc1","Myh6","Ttn","Acta2","Myl9","Tnnt2","Myl4",
              "Tpm1","Myl2","Acta1","Actn2","Atp2a2","Bmp2","Cdh2","Csrp2",
              "Csrp3","Mef2c","Mybpc3","Nkx2.5","Qki","Slc8a1","Tbx2","Tmod1",
              "Tpm1","Ldb3","Hacd1","Pdlim5","Cacna2d2","Popdc2","Asb2","Smarcd3",
              "Nebl","Myo18b","Lmod1","Myh7","Myocd","Sorbs2","Cnn1","Cnn2","Myo18a",
              "Myl6","Ryr2","Tpm2","Dstn","Pdlim7","Palld","Svil","Mical2","Fitm1",
              "Ank3","Cck","Afdn","Nexn","Ccdc141","Ehd4", ##心肌
              "Ptn","Dlk1","Peg3","Meg3", "Hmgb2","Birc5","Cdk1","Cdk4","Ezh2","Gas1","Hmgb1",
              "Hmga2","Lmnb1","Mki67","Meis2","Rrm2","Cks1b","Dact1","Cks2","Anp32b","Ube2c",
              "Cenpf","Smc2","Tgfb2","Wnt4","Smc4","Stmn1","Ran","Tubb5","Nasp","Usp29","Hoxb4",
              "Top2a","Gpc3")
seurat.object_copy3<-ScaleData(seurat.object_copy2, features = cg_markers)
DoHeatmap(seurat.object_copy3,
          cg_markers,
          size=3,
          group.by = "highORlow") 
DoHeatmap(seurat.object_copy3,
          cg_markers,　
          group.by = "highORlow",
          label=F)+
  scale_fill_gradientn(colors = c("lightblue","white","firebrick3"))
#scale_fill_viridis()


######差异基因
markers <- FindMarkers(seurat.object_copy2, ident.1 = "double_pos", 
                       ident.2 = "Hcn4_pos",
                       group.by = 'highORlow',
                       logfc.threshold = 0.25)
head(x = markers)
write.csv(markers,"d13-double_pos_vs_shox2_pos.csv")

##筛选后作图
cg_markers=markers[abs(markers$avg_log2FC) >0.5,]
dim(cg_markers)
##这里会遇到报错，如果基因非高变基因的话
##故尝试myda2<- ScaleData(mydata, features = all.genes)
seurat.object_copy3<-ScaleData(seurat.object_copy2, features = rownames(cg_markers))
DoHeatmap(seurat.object_copy3,
          rownames(cg_markers) ,
          size=3,
          group.by = "highORlow") 

DoHeatmap(seurat.object_copy3,
          rownames(cg_markers),　
          group.by = "highORlow",
          label=F)+
  scale_fill_gradientn(colors = c("lightblue","white","firebrick3"))
#scale_fill_viridis()



##shox2+ & hcn4+
select_cells <- WhichCells(seurat.object_copy, slot = 'counts', expression = Hcn4 > 0 & Tbx3 > 0 & Shox2 > 0 ) 
select_obj <- subset(seurat.object_copy, cells = select_cells)


SAN_head=colnames(subset(x = select_obj, subset = `Tbx18` > 0 & `Shox2` > 0 & `Nkx2-5` <= 0 & Tbx3 >0, slot = 'counts'))
SAN_tail =colnames(subset(x = select_obj, subset = Tbx18 <= 0 & Shox2 > 0 & `Nkx2-5` > 0 & Tbx3 >0, slot = 'counts'))
SAN_TZ=colnames(subset(x = select_obj, subset = Shox2 <= 0 & Nppa > 0 & `Nkx2-5` > 0 & Tbx3 >0, slot = 'counts'))

highORlow=ifelse(colnames(select_obj) %in% SAN_head,'SAN_head',
                 ifelse(colnames(select_obj) %in% SAN_tail,'SAN_tail',
                        ifelse(colnames(select_obj) %in% SAN_TZ,'SAN_TZ','others')))
table(highORlow)
select_obj@meta.data$highORlow=highORlow
DimPlot(select_obj, reduction = "umap", label = T, group.by = "highORlow", pt.size = 0.5)
DimPlot(select_obj, reduction = "umap", label = T, split.by = "highORlow",pt.size = 0.5)























DimPlot(select_obj, reduction = "umap", label=T) 

markers.to.plot <- c("Tbx3","Tbx18","Shox2","Nkx2-5","Nppa")
FeaturePlot(select_obj, features = markers.to.plot,label = TRUE)
FeaturePlot(select_obj, features = markers.to.plot,split.by = "sample",label = F)
p<-FeaturePlot(scRNA_harmony, features = c("Hcn4","Shox2"), split.by = "sample",label = TRUE)
p+scale_color_viridis()
markers.to.plot <- c("Hcn4","Shox2")
##cividis,magma,inferno,viridis(默认)
DotPlot(select_obj, features = rev(markers.to.plot), dot.scale = 8) + RotatedAxis()+scale_color_viridis(option="cividis")
DotPlot(select_obj, features = rev(markers.to.plot), dot.scale = 8) + scale_color_viridis()
DotPlot(select_obj, features = rev(markers.to.plot)) + scale_color_viridis()

VlnPlot(select_obj, features = markers.to.plot,log = F,slot = "data",pt.size = 0)
VlnPlot(select_obj, features = markers.to.plot,log = F,slot = "data",pt.size = 0.1)
##尝试对细胞进行分类并命名

