

setwd("E:/workspace/keyan/GC_tme/Article_results/7.viper_activity_kim2018/")

library("viper")

eset = read.table("ICB.Kim2018_Pembrolizumab_Gastric.RNAseq.txt",
                  sep="\t",header = TRUE,check.names = FALSE,row.names = 1)

regulons <- aracne2regulon(afile="diff_network_sub.txt",eset = as.matrix(eset),format="3col")
save(regulons,file="regul.RData")

vpres <- viper(eset, regulons, verbose = FALSE)

save(vpres,file="viper.RData")


result = as.data.frame(vpres)
write.table(cbind(gene=rownames(result),result),file="ICB_kim2018_viper_activity_proteins.txt",
            sep="\t",row.names = FALSE,quote = FALSE)


