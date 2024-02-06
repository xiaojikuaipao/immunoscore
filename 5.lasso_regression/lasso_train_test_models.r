#!/usr/bin/Rscript


setwd("E:/workspace/keyan/GC_tme/7.TCGA-STAD_immune_cal/4.5.arace_ap_network_viper/diff_exp_sub_network/diff_cal_protein_activity/Venn_Lasso")
suppressMessages(library('glmnet'))
library("tidyverse")
library("pROC")



expf = "viper_activity_proteins.txt"
grpf = "immune_cluster_ordered.txt"


family = "multinomial"


modeldata = try(read.table(expf,sep="\t",row.names=1, header=TRUE,check.names=F,encoding="utf-8"))

groupdata = try(read.table(grpf,sep="\t",row.names =1,header=TRUE,check.names=F,encoding="utf-8"))

trainx = t(modeldata[,rownames(groupdata)])
trainy = groupdata$cluster


fitcv <- cv.glmnet(trainx,trainy,family=family,nfold=10,intercept = F, alpha=1)

png("LASSO_coefficient_CV.png",width=4000,height=2000,res=300)
plot(fitcv)
dev.off()

png("LASSO_coefficient.png",width=4000,height=2000,res=300)
plot(fitcv$glmnet.fit,xvar="lambda",cex=1.5,cex.axis=1.5,cex.lab=1.5)
dev.off()

pdf("LASSO_coefficient_CV.pdf")
plot(fitcv)
dev.off()

pdf("LASSO_coefficient.pdf")
plot(fitcv$glmnet.fit,xvar="lambda",cex=1.5,cex.axis=1.5,cex.lab=1.5)
dev.off()

print("model---")
print(fitcv$lambda.min,fitcv$lambda.1se)
write.table(file="link.txt",data.frame(min=fitcv$lambda.min,se=fitcv$lambda.1se),sep="\t")
	
fitcv_coef <- coef(fitcv$glmnet.fit,s=fitcv$lambda.min)

result = NA
for( i in names(fitcv_coef)){
  e = as.data.frame(as.matrix(fitcv_coef[[i]]))
  result = cbind(e,result)
}

nresult = result[-1,-5]
colnames(nresult) = c("cluster1","cluster2","cluster3","cluster4")
resultf = 'LASSO_coefficient_list.xls'

write.table(cbind(gene=rownames(nresult),nresult), resultf,sep='\t',quote=F,row.names = FALSE)
  
 	
filtered_result = nresult[apply(nresult,1,sum) >0,] 
resultf_filter = 'LASSO_coefficient_list_filter.xls'
write.table(cbind(gene=rownames(filtered_result),filtered_result), resultf_filter,sep='\t',quote=F,row.names = FALSE)


filtered_modelf = "viper_activity_proteins_filter_lasso_keygene.txt"
modeldata_filtered = modeldata[rownames(filtered_result),]
write.table(cbind(gene=rownames(modeldata_filtered),modeldata_filtered), filtered_modelf,sep='\t',quote=F,row.names = FALSE)


