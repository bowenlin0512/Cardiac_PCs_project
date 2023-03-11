###monocle analysis
library(monocle)
library(ggsci)
library(ggpubr) 
library(patchwork)
##prepare the input
SAN_sub = seurat.object_copy[,seurat.object_copy@meta.data$cell_type %in% c("Nascent mesoderm",
                                                                            "Cardiac progenitor cell","SHF","Pacemaker progenitor cell",
                                                                            "Immature cardiomyocyte","Cardiac fibroblast","SAN")]
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(SAN_sub$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = SAN_sub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#cds <- detectGenes(cds, min_expr = 0.5) 
cds <- detectGenes(cds, min_expr = 0.5) #add num_cells_expressed 
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 100))

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type",cores= 1) 
head(diff) 
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),] 
head(deg)
##save the train.monocle.DEG.xls
write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
deg<-read.csv("train.monocle.DEG.csv",header = T)
ordergene <- deg$gene_short_name 
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene) 
pdf("train.ordergenes.pdf") 
plot_ordering_genes(cds)

##Appliying the top genes
#ordergene <- row.names(deg)[order(deg$qval)][1:2000]
ordergene <- deg$gene_short_name[order(deg$qval)][1:2000]
##save the cds file
#saveRDS(cds, file = "mococle_before.rds")

##Step 2. 降维
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
## Trajectory step 2: reduce data dimensionality----
# cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
## Trajectory step 3: order cells along the trajectory ----
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 5) 
#saveRDS(cds, file = "mococle_final.rds")
cds<-readRDS("mococle.rds")

##plots
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "cell_type")
plot_cell_trajectory(cds, color_by = "State")
# "State" is just Monocle's term for the segment of the tree.
plot_cell_trajectory(cds, color_by="State") +
  facet_wrap(~State, nrow=1)
cds<- orderCells(cds, root_state = 1)
plot_cell_trajectory(cds, color_by = "Pseudotime")
p1<-plot_cell_trajectory(cds, color_by = "State",cell_size = 0.75)+scale_color_nejm(alpha = 1) 
plot_cell_trajectory(cds, color_by = "State",cell_size = 0.75)+facet_wrap(~State) 
plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 0.75) 
plot_cell_trajectory(cds, color_by = "seurat_clusters",cell_size = 0.75) +scale_color_manual(values=col) 
plot_cell_trajectory(cds, color_by = "seurat_clusters", cell_size = 0.75) + facet_wrap (~seurat_clusters, nrow =Cluster_num) 
plot_cell_trajectory(cds, color_by = "sample",cell_size = 0.75) 
p3<-plot_cell_trajectory(cds,color_by="cell_type", size=1,show_backbone=TRUE)+scale_color_lancet(alpha = 1)
p2<-plot(plot_cell_trajectory(cds, show_cell_names = F, color_by = "Pseudotime") + scale_color_viridis_c())
##
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "cell_type") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend 
  scale_color_manual(values = colour) 
p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "cell_type")+ 
  scale_color_manual(values = colour) + theme(legend.title = element_blank()) 
p1|p2|p3
p2|p3|p1
##Density

df <- pData(cds) 
View(df) 
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) + 
  geom_density(bw=0.5,size=1,alpha = 0.5)+scale_color_npg()

##set the genes
s.genes <- c("Zfp503","Lbh","Bmp2","Socs2") 
s.genes <- c("Zfp503")
s.genes <- c("Tbx18")
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State") 
p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State") 
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")+scale_color_npg()
plotc <- p1|p2|p3 
ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 8)
pData(cds)$Shox2 = log2( exprs(cds)['Shox2',]+1)
plot_cell_trajectory(cds, color_by = "Shox2")  + scale_color_viridis()

###show the expreesion of genes
colnames(pData(cds))
pData(cds)$Gsc = log2( exprs(cds)['Gsc',]+1)
p1=plot_cell_trajectory(cds, color_by = "Gsc",cell_size = 0.5)  + scale_color_viridis()
pData(cds)$Hcn4 = log2( exprs(cds)['Hcn4',]+1)
p2=plot_cell_trajectory(cds, color_by = "Hcn4",cell_size = 0.5)  + scale_color_viridis()
pData(cds)$Shox2 = log2(exprs(cds)['Shox2',]+1)
p3=plot_cell_trajectory(cds, color_by = "Shox2",cell_size = 0.5)+ scale_color_viridis()
pData(cds)$Col1a1 = log2( exprs(cds)['Col1a1',]+1)
p4=plot_cell_trajectory(cds, color_by = "Col1a1",cell_size = 0.5)  + scale_color_viridis()
pData(cds)$Tbx18 = log2( exprs(cds)['Tbx18',]+1)
p5=plot_cell_trajectory(cds, color_by = "Tbx18",cell_size = 0.5)  + scale_color_viridis()
pData(cds)$Isl1 = log2( exprs(cds)['Isl1',]+1)
p6=plot_cell_trajectory(cds, color_by = "Isl1",cell_size = 0.5)  + scale_color_viridis()
p1+p2+p3+p5+p6+p4+plot_layout(nrow = 2)

##for Zfp503
pData(cds)$Zfp503 = log2( exprs(cds)['Zfp503',]+1)
plot_cell_trajectory(cds, color_by = "Zfp503",cell_size = 0.5)  + scale_color_viridis()

##heatmap
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
cds_sub <- cds[disp.genes,]
plot_cell_trajectory(cds_sub, color_by = "State")
plot_cell_trajectory(cds_sub, color_by = "cell_type")
plot_cell_trajectory(cds_sub, color_by = "Pseudotime",cell_size = 0.75) ##根据拟时间值着色
beam_res <- BEAM(cds_sub, branch_point = 1, cores = 1)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
cds_sub_beam <- cds_sub[row.names(subset(beam_res, qval < 1e-4)),]
pdf("genes_branched_heatmap.pdf",width = 9,height = 5)
plot_genes_branched_heatmap(cds_sub_beam,  branch_point = 1, num_clusters = 5, 
                            hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                            branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                            use_gene_short_name = T,
                            show_rownames = F,
                            return_heatmap = T #是否返回一些重要信息
)
dev.off()
p<-plot_genes_branched_heatmap(cds_sub_beam,  branch_point = 1, num_clusters = 5, 
                               hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                               branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                               use_gene_short_name = T,
                               show_rownames = F,
                               return_heatmap = T #是否返回一些重要信息
)
pdf("genes_branched_heatmap.pdf",width = 9,height = 5)
plot_genes_branched_heatmap(cds_sub_beam, branch_point = 1,
                            branch_states = NULL, branch_labels = c("Cell fate 1", "Cell fate 2"),
                            cluster_rows = TRUE, hclust_method = "ward.D2", num_clusters = 6,
                            hmcols = NULL, branch_colors = c("#979797", "#F05662", "#7990C8"),
                            add_annotation_row = NULL, add_annotation_col = NULL,
                            show_rownames = FALSE, use_gene_short_name = TRUE, scale_max = 3,
                            scale_min = -3, norm_method = c("log", "vstExprs"),
                            return_heatmap = F, cores = 1)
p<-plot_genes_branched_heatmap(cds_sub_beam, branch_point = 1,
                               branch_states = NULL, branch_labels = c("Cell fate 1", "Cell fate 2"),
                               cluster_rows = TRUE, hclust_method = "ward.D2", num_clusters = 6,
                               hmcols = NULL, branch_colors = c("#979797", "#F05662", "#7990C8"),
                               add_annotation_row = NULL, add_annotation_col = NULL,
                               show_rownames = FALSE, use_gene_short_name = TRUE, scale_max = 3,
                               scale_min = -3, norm_method = c("log", "vstExprs"),
                               return_heatmap = T, cores = 1)

dev.off()
saveRDS(mycds, file="mycds.rds")
test_genes=c("Ets1","Rspo3","Ptch1","Robo2")
test_genes=c("Otx2","Zfp503","Nkx2-5","Shox2","Klf6","Cenpa","Nfat5","Tbx20","Foxp1","Zeb1")
test_genes=c("Zfp503")
pdf("genes_branched_pseudotime3.pdf",width = 9,height = 5)
plot_genes_branched_pseudotime(cds_sub_beam[test_genes,],
                               branch_point = 1,
                               color_by = "cell_type",
                               cell_size=2,
                               ncol = 4)+scale_color_lancet(alpha = 1)
dev.off()