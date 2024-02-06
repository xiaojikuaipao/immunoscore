
rm(list=ls())

library("ggplot2")
library("patchwork")
library("dplyr")
library("tidyverse")
library("openxlsx")
library("rstatix")

setwd("E:/workspace/keyan/GC_tme/Article_results/8.ICB_kim2018_roc/")

library("survival")
library("survminer")

exp = read.table("ICB_kim2018_viper_activity_proteins.txt",sep="\t",
                 header = TRUE,row.names = 1,
                 check.names = FALSE)
genes = c("MYB","NFE2L3","SETBP1","TFAP2A")
nexp = exp[genes,]
pc <- prcomp(t(nexp))
value <- pc$x[, 1] * sign(cor(pc$x[, 1], colMeans(nexp,na.rm = T)))

result = as.data.frame(value)

result_final = data.frame(sample=rownames(result),result)
write.table(result_final,
            file="sample_PCA_score.txt",row.names = FALSE,sep="\t",
            quote = FALSE)

phe = read.table("ICB.Kim2018_Pembrolizumab_Gastric.clinical.txt",sep="\t",header = TRUE)
rownames(phe) = phe$patient

nphe = phe[colnames(nexp),]
#q1 <- quantile(nresult$value, 0.25)
#median_value <- quantile(nresult$value, 0.5)
#q3 <- quantile(nresult$value, 0.75)


result$group = nphe$Response

result$group = as.factor(result$group)
result$value = as.numeric(result$value)

result = as.data.frame(result)

stat.test <- result %>% 
  t_test(value~group)

stat.test <- stat.test %>% 
  add_xy_position(x='group',dodge = 1)
#add_y_position()

p = ggboxplot(result, x = "group", y = "value", color = "group", 
          notch = FALSE, outlier.shape = 1, outlier.size = 0.5, 
          add = "jitter",  # 添加 jitter 以更好地显示数据分布
          add.params = list(alpha = 1)) +  # 设置 jitter 的透明度
  #stat_compare_means(aes(group=group),method = "wilcox.test",label = "p.format")+
  stat_pvalue_manual( stat.test,  label = "p") +
  theme(
    axis.title.x = element_blank(),  # 去掉 x 轴标题
    axis.text.x = element_text(size = 15, face = "bold",angle=45, hjust=1),
    panel.background = element_blank()  # 去掉图布背景
  )

pdf("Kim_validate_score.pdf",width = 4,height = 7)
p
dev.off()

png("Kim_validate_score.png",width=4,height=7,units="in",res=300,pointsize=18)
p
dev.off()


write.table(cbind(gene=rownames(result),result),file="sample_PCA_score_group.txt",sep="\t",row.names = FALSE,quote = FALSE)

############

library(pROC)

pdf("ROC.pdf",width = 6,height = 6)
rocobj = roc(result$group,as.numeric(result$value))
plot(rocobj,print.auc=TRUE,auc.polygon=FALSE,grid=c(0.1,0.2),col="red",
     print.thres=TRUE,main="")
dev.off()
pars = coords(rocobj, "best", ret=c("threshold", "specificity", "sensitivity"))
vs = power.roc.test(rocobj)
sig = vs$sig.level
cis = as.character(ci.auc(rocobj))
upci = cis[1]
downci = cis[3]
aucs = cis[2]

outdata = data.frame("PCA_score",threshold=as.character(pars[1]),
                     specificity=as.character(pars[2]),
                     sensitivity=as.character(pars[3]),
                     auc=aucs,upci=upci,downci=downci)
write.table(outdata,file="ROC_static.txt")









