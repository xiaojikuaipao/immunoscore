


library("pheatmap")

setwd("E:/workspace/keyan/GC_tme/Article_results/3.Limma_screen_biomarker/Venn_heatmap/")

exp = read.table("TCGA-STAD_fpkm_cluster_diff.txt",sep="\t",header = TRUE,check.names = FALSE,row.names = 1)
phe = read.table("immune_cluster_ordered.txt",sep="\t",header = TRUE,row.names =1,check.names = FALSE)

nexp = exp[,rownames(phe)]
colnames(phe) = "Clusters"


p = pheatmap(nexp,annotation_col = phe,cluster_rows = TRUE,cluster_cols = FALSE,
         color = colorRampPalette(c("navy","yellow","firebrick3"))(10),
         breaks=seq(0,10,1),
         show_rownames = FALSE,show_colnames = FALSE)

pdf("cluster_heatmap.pdf")
print(p)
dev.off()

png("cluster_heatmap.png")
p

dev.off()

