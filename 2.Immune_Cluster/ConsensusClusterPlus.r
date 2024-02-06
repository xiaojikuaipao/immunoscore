
setwd("E:/workspace/keyan/GC_tme/Article_results/2.Immune_Cluster/")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#  BiocManager::install("ConsensusClusterPlus")

library(RColorBrewer)
library(pheatmap)
library("ConsensusClusterPlus")

immune=read.table("CIBERSORT-Results.txt",sep = "\t",
                  header = TRUE,check.names = FALSE,row.names = 1)
cols = ncol(immune)
immune_final = immune[,-c(cols,cols-1,cols-2)]
immune_ddist = t(immune_final)

ConsensusClusterObject <- ConsensusClusterPlus(immune_ddist,                                                                
                                               maxK = 6,                # set K
                                               pItem = 0.8,
                                               reps=5000,              # set repeats
                                               clusterAlg = "hc",       # clustering Algorithm : Hierarchiocal clustering
                                               innerLinkage = "ward.D2",  # for color coding the clusters use tmyPal = ...
                                               finalLinkage = "complete",
                                               plot = 'pdf',      # write resut to pdf (Alt:png)
                                               writeTable = TRUE,
                                               title="immune_cluster",
                                               verbose = TRUE)

save(ConsensusClusterObject,file="ConsensusClusterObject.RData")

optimal_k = "4"
classf = paste0("immune_cluster/","immune.k=",optimal_k,".consensusClass.csv")
class_df = read.table(classf,sep=",",header = FALSE)
rownames(class_df) = class_df$V1

outclass = data.frame(sample=class_df$V1,cluster=paste0("Cluster",class_df$V2))
write.table(outclass,file="consensus_cluster.txt",sep="\t",row.names = FALSE,quote = FALSE)

outclass_sorted = outclass[order(outclass$cluster),]
immune_ddist_sorted = immune_ddist[,outclass_sorted$sample]
cluster_ano = data.frame(Clusters=outclass_sorted$cluster,row.names = outclass_sorted$sample)

write.table(cbind(cell=rownames(immune_ddist_sorted),immune_ddist_sorted),"immune_ordered_cell_percent.txt",sep="\t",row.names = FALSE,quote = FALSE)
write.table(outclass_sorted,"immune_cluster_ordered.txt",sep="\t",row.names = FALSE,quote = FALSE)

gap_cols = c()
e = 0
for( i in paste0("Cluster",as.character(seq(1,4,1)))){
  v = nrow(outclass_sorted[outclass_sorted$cluster == i,])
  print(v)
  gap_cols = append(gap_cols,e+v)
  e = e+v
}


p = pheatmap(immune_ddist_sorted,
             cluster_cols = F,
             cluster_rows = F,
             annotation_col = cluster_ano,
             show_colnames = F,
             show_rownames = T,
             color = colorRampPalette(c("navy","white","firebrick3"))(100),
             #breaks = c(0,10),
             fontsize = 16,
             scale = "column",
             #breaks = seq(-4,4,1)
             gaps_col = gap_cols,
             cellwidth = 0.8           # 调整单元格宽
             
)

pdf("immune_cluster.pdf",width = 10,height = 10)
print(p)
dev.off()

png("immune_cluster.png",width = 1000,height = 1000)
p
dev.off()



