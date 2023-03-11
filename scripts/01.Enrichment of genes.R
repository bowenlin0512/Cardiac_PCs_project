##Enrichment (GO or KEGG) of genes
##GO
gene_group=p$annotation_row
gene_group$gene=rownames(gene_group)
c <- bitr(gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
d<-merge(gene_group,c,by.x="gene",by.y="SYMBOL")
write.csv(d,"branched_pseudotime_genes.csv")
library(clusterProfiler)
library(org.Mm.eg.db)
#####compare kegg
kegg <- compareCluster(ENTREZID~Cluster, 
                       data=d,
                       fun="enrichKEGG", 
                       organism="mmu", 
                       #organism="rno"
                       pvalueCutoff=0.5,
                       pAdjustMethod = "BH", 
                       qvalueCutoff = 0.5)
dotplot(kegg,showCategory=15,includeAll=TRUE)
y2 <- setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(y2,"beam_kegg.csv")
#####compare GO
c<-na.omit(c)
go <- compareCluster(ENTREZID~Cluster, 
                     data=d,
                     fun="enrichGO", 
                     ont= "BP",
                     OrgDb="org.Mm.eg.db", 
                     #OrgDb="org.Rn.eg.db",
                     pvalueCutoff=0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1)
ego2 <- simplify(go,cutoff=0.7,by="p.adjust",select_fun=min)
dotplot(ego2,showCategory=5,includeAll=TRUE)

allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])
write.csv(allcluster_go,"beam_GO.csv")



test_genes=c("TFF3","GUCA2B")
pdf("genes_branched_pseudotime.pdf",width = 9,height = 4)
plot_genes_branched_pseudotime(test[test_genes,],
                               branch_point = 1,
                               color_by = "celltype",
                               cell_size=2,
                               ncol = 2)
dev.off()