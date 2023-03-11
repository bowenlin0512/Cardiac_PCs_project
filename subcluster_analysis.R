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
select_cells <- WhichCells(seurat.object_copy, slot = 'counts', expression = Tbx18 > 0) 
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
