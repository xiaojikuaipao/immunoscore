

setwd("E:/workspace/keyan/GC_tme/Article_results/4.aracne_viper_activity/")

library("viper")

eset = read.table("TCGA-STAD_fpkm_cluster_diff.txt",
                  sep="\t",header = TRUE,check.names = FALSE,row.names = 1)

regulons <- aracne2regulon(afile="diff_network_sub.txt",eset = as.matrix(eset),format="3col")
save(regulons,file="regul.RData")

vpres <- viper(eset, regulons, verbose = FALSE)

save(vpres,file="viper.RData")


result = as.data.frame(vpres)
write.table(cbind(gene=rownames(result),result),file="viper_activity_proteins.txt",sep="\t",row.names = FALSE,quote = FALSE)


