##scRNA SAN
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(RColorBrewer)
library(viridis)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
#BiocManager::install("monocle")
#remotes::install_github("satijalab/seurat")
library(DoubletFinder)
library(monocle)
setwd("/home/bowenlin/disk1/project/yang/lzy/SAN_scRNA-seq/")
scRNA_harmony<-readRDS("scRNA_harmony_before_definaiton.rds")
DimPlot(scRNA_harmony, reduction = "umap", label = T, pt.size = 0.5)+scale_color_simpsons()

##Dimplot change clolors
#--------------------------------------------------
library(RColorBrewer)
library(paletteer)
library(ggsci)
my_cols <- c('3'='#F68282','15'='#31C53F','5'='#1FA195','1'='#B95FBB','13'='#D4D915',
             '14'='#28CECA','9'='#ff9a36','8'='#2FF18B','11'='#aeadb3','6'='#faf4cf',
             '2'='#CCB1F1','12'='#25aff5','7'='#A4DFF2','4'='#4B4BF7','16'='#AC8F14',
             '10'='#E6C122')
#选择Set1的前5个和Set2的后4个使用。
pal1<-brewer.pal(9,'Set1')
pal2<-brewer.pal(8,'Set2')
pal3<-brewer.pal(12,'Set3')
mycolor<-c(pal1[1:9],pal2[1:8],pal3[1:6])
DimPlot(scRNA_harmony, label = T, cols = cell_type_cols, #设置颜色 
        pt.size = 1, #设置点的大小 
        repel = T) #标注有点挤，repel=T可以让排列更加合理 
#--------------------------------------------------      
        


dir.create("pseudotime")
combined.6 <- scRNA_harmony
combined.7 <- RenameIdents(combined.6, 
                           '0'='E4.75-5',
                           '1'='E4.75-5',
                           '2'='2',
                           '3'='3',
                           '4'='4',
                           '5'='5',
                           '6'='6',
                           '7'='7',
                           '8'='core',
                           '9'='Definitive Endoderm',
                           '10'='10',
                           '11'='11',
                           '12'='12',
                           '13'='13',
                           '14'='14',
                           '15'='15',
                           '16'='16',
                           '17'='17',
                           '18'='18',
                           '19'='posterior SHF',
                           '20'='20',
                           '21'='Hepatic Endoderm',
                           '22'='Endothelial cells');
levels(combined.7)
DimPlot(combined.7, reduction = "umap", label = T, pt.size = 0.8)
VlnPlot(combined.7, features = select_genes, pt.size=0, ncol=3)
DotPlot(combined.7, features = select_genes,cols = c("grey", "red"),) + coord_flip()
DotPlot(combined.7, features = all_markers,cols = c("grey", "red"),) 
dev.off()
save(combined.7,file="RunUMAP_rename.Robj") 

###########################monocle
dir.create("pseudotime")
combined.7@meta.data$celltype<-combined.7@active.ident
Cells.sub <- subset(combined.7@meta.data, celltype=="10"|celltype=="core"|celltype=="13")
SAN.mergesub <- subset(scRNA_harmony, cells=row.names(Cells.sub))

data <- as(as.matrix(SAN.mergesub@assays$RNA@counts), 'sparseMatrix')
# count矩阵
pd <- new('AnnotatedDataFrame', data = SAN.mergesub@meta.data)
# meta表转成特定格式
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# 基因名表转成特定格式
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
#expressionFamily参数用于指定表达矩阵的数据类型，有几个选项可以选择：
#稀疏矩阵用negbinomial.size()，FPKM值用tobit()，logFPKM值用gaussianff()
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

##使用clusters差异表达基因
#diff.genes <- read.csv('subcluster/diff_genes_wilcox.csv')
#diff.genes <- subset(diff.genes,p_val_adj<0.01)$gene
#mycds <- setOrderingFilter(mycds, diff.genes)
#p1 <- plot_ordering_genes(mycds)
##使用seurat选择的高变基因
#var.genes <- VariableFeatures(SAN.mergesub)
#mycds <- setOrderingFilter(mycds, var.genes)
#p2 <- plot_ordering_genes(mycds)
##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
##结果对比
#p1|p2|p3
p3

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("pseudotime/State.pdf", plot = plot1, width = 6, height = 5)
ggsave("pseudotime/State.png", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("pseudotime/Cluster.pdf", plot = plot2, width = 6, height = 5)
ggsave("pseudotime/Cluster.png", plot = plot2, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("pseudotime/Pseudotime.pdf", plot = plot3, width = 6, height = 5)
ggsave("pseudotime/Pseudotime.png", plot = plot3, width = 6, height = 5)
##合并作图
plotc <- plot1|plot2|plot3
ggsave("pseudotime/Combination.pdf", plot = plotc, width = 10, height = 3.5)
ggsave("pseudotime/Combination.png", plot = plotc, width = 10, height = 3.5)
##保存结果
write.csv(pData(mycds), "pseudotime/pseudotime.csv")

##轨迹分面展示
p1 <- plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~State, nrow = 1)
p2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters, nrow = 1)
plotc <- p1/p2
ggsave("pseudotime/trajectory_facet.png", plot = plotc, width = 6, height = 5)


s.genes <- c("Shox2","Hcn4","Nanog")
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("pseudotime/genes_visual.png", plot = plotc, width = 8, height = 4.5)

#cluster差异基因
#diff.genes <- read.csv('subcluster/diff_genes_wilcox.csv')
#sig_diff.genes <- subset(diff.genes,p_val_adj<0.0001&abs(avg_logFC)>0.75)$gene
#sig_diff.genes <- unique(as.character(sig_diff.genes))
#diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
#                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
#sig_gene_names <- row.names(subset(diff_test, qval < 0.01))
#p1 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=3,
#                             show_rownames=T, return_heatmap=T)
#ggsave("pseudotime/pseudotime_heatmap1.png", plot = p1, width = 5, height = 8)
#高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))
p2 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime/pseudotime_heatmap2.png", plot = p2, width = 5, height = 10)
#高变基因2
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 2&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))
p2 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime/pseudotime_heatmap2.png", plot = p2, width = 5, height = 10)
####BEAM
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- mycds[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 3, show_rownames = T)
saveRDS(mycds, file="mycds.rds")



#根据Shox2及Hcn4表达量进行差异分析
#提取d13
Cells.sub <- subset(scRNA_harmony@meta.data, sample=="d13")
SAN_d13 <- subset(scRNA_harmony, cells=row.names(Cells.sub))
##查看提取结果
DimPlot(SAN_d13, reduction = "umap", label = T, pt.size = 0.5)

##提取d13 Shox2及Hcn4 细胞
double_pos=colnames(subset(x = SAN_d13, subset = Shox2 > 0  & Hcn4 > 0 , slot = 'counts'))
shox2_pos=colnames(subset(x = SAN_d13, subset = Shox2 > 0 & Hcn4 <= 0 , slot = 'counts'))
Hcn4_pos=colnames(subset(x = SAN_d13, subset = Shox2 <= 0 & Hcn4 > 0 , slot = 'counts'))

highORlow=ifelse(colnames(SAN_d13) %in% double_pos,'double_pos',
                 ifelse(colnames(SAN_d13) %in% shox2_pos,'shox2_pos',
                 ifelse(colnames(SAN_d13) %in% Hcn4_pos,'Hcn4_pos','double_neg')))
table(highORlow)
SAN_d13@meta.data$highORlow=highORlow
markers <- FindMarkers(SAN_d13, ident.1 = "double_pos", 
                       ident.2 = "shox2_pos",
                       group.by = 'highORlow',
                       logfc.threshold = 0.25)
head(x = markers)
write.csv(markers,"d13-double_pos_vs_shox2_pos.csv")
##筛选后作图
cg_markers=markers[abs(markers$avg_log2FC) >0.5,]
dim(cg_markers)
DoHeatmap(SAN_d13,
          rownames(cg_markers) ,
          size=3,
          group.by = "highORlow") 

DoHeatmap(SAN_d13,
          rownames(cg_markers),　
          group.by = "highORlow",
          label=F)+
  scale_fill_gradientn(colors = c("lightblue","white","firebrick3"))
   #scale_fill_viridis()

##类似地，提取d8的细胞
#----------------------
Cells.sub <- subset(scRNA_harmony@meta.data, sample=="d8")
SAN_d8 <- subset(scRNA_harmony, cells=row.names(Cells.sub))
##查看提取结果
DimPlot(SAN_d8, reduction = "umap", label = T, pt.size = 0.5)

##提取d8 Shox2及Hcn4 细胞
double_pos=colnames(subset(x = SAN_d8, subset = Shox2 > 0  & Hcn4 > 0 , slot = 'counts'))
shox2_pos=colnames(subset(x = SAN_d8, subset = Shox2 > 0 & Hcn4 <= 0 , slot = 'counts'))
Hcn4_pos=colnames(subset(x = SAN_d8, subset = Shox2 <= 0 & Hcn4 > 0 , slot = 'counts'))

highORlow=ifelse(colnames(SAN_d8) %in% double_pos,'double_pos',
                 ifelse(colnames(SAN_d8) %in% shox2_pos,'shox2_pos',
                        ifelse(colnames(SAN_d8) %in% Hcn4_pos,'Hcn4_pos','double_neg')))
table(highORlow)
SAN_d8@meta.data$highORlow=highORlow
markers <- FindMarkers(SAN_d8, ident.1 = "double_pos", 
                       ident.2 = "Hcn4_pos",
                       group.by = 'highORlow',
                       logfc.threshold = 0.25)
head(x = markers)
write.csv(markers,"d8-double_pos_vs_Hcn4_pos.csv")
##筛选后作图
cg_markers=markers[abs(markers$avg_log2FC) >0.4,]
dim(cg_markers)
DoHeatmap(SAN_d8,
          rownames(cg_markers) ,
          size=3,
          group.by = "highORlow") 

DoHeatmap(SAN_d8,
          rownames(cg_markers),　
          group.by = "highORlow")+
  scale_fill_gradientn(colors = c("lightblue","white","firebrick3"))
  #scale_fill_viridis()




