

rm(list=ls())
setwd("E:/workspace/keyan/GC_tme/Article_results/6.survival_plot_stat/")

library("survival")
library("survminer")

exp = read.table("viper_activity_proteins_filter_lasso_keygene.txt",sep="\t",
                 header = TRUE,row.names = 1,
                 check.names = FALSE)


surv = read.table("TCGA-STAD.survival.tsv",sep="\t",header = TRUE,
                  row.names = 1,check.names = FALSE)

common_name = intersect(colnames(exp),rownames(surv))

nsurv = surv[common_name,]
nexp = exp[,common_name ]
colnames(nsurv) = c("statu","patient","time")

if(!file.exists("survivalplot")){
  dir.create("survivalplot")
}
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5,size =20,face="bold")
    )
}

out = data.frame(gene=NA,pvalue=NA)
for(g in rownames(nexp)){
  print(nexp[g,])
  median_value = median(as.numeric(nexp[g,]))
  nsurv$group = ifelse(as.numeric(nexp[g,])>median_value,"High","Lower")
  fit <- survfit(Surv(time,statu) ~ group,  data = nsurv) 
  
 
  pdf(paste0("survivalplot/",g,"_survival_plot.pdf"))
  myplot = ggsurvplot(fit, # 创建的拟合对象
             title = g,
             data = nsurv,  # 指定变量数据来源"
             conf.int = TRUE, # 显示置信区间
             xlab="Time(days)",
             ylab = "Overall survival",
             pval = TRUE, # 添加P值
             pval.method = TRUE,
             risk.table = FALSE, # 绘制累计风险曲线
             surv.median.line = "hv", # 添加中位生存时间线
             add.all = FALSE, # 添加总患者生存曲线
             palette = "nejm",
             ggtheme=custom_theme(),
             legend = c(0.8,0.75),
             legend.labs=c("High","Low"),
             legend.title="Activity")
    
  
  print(myplot,newpage = FALSE)
  dev.off()
  #ggsave(paste0("survivalplot/",g,"_survival_plot.pdf"),print(myplot), device = "pdf") 
  surv_diff_result <- survdiff(Surv(time,statu) ~ group,  data = nsurv)
  p_value_logrank <- 1 - pchisq(surv_diff_result$chisq, df = 1)
  out = rbind(out,data.frame(gene=g,pvalue=p_value_logrank))

}
write.table(out[-1,],file="gene_survival_pvalue_stats.txt",sep="\t",quote = FALSE,row.names = FALSE)

pvalue_stat = out[-1,]
pvalue_sig_df = pvalue_stat[pvalue_stat$pvalue < 0.05,]
pvalue_expf = nexp[pvalue_sig_df$gene,]
write.table(cbind(gene=rownames(pvalue_expf),pvalue_expf),file="viper_activity_proteins_filter_lasso_keygene_filter_survival.txt",
            sep="\t",row.names = FALSE,quote = FALSE)

#for(g in rownames(pvalue_expf)){
#  median_value = median(as.numeric(pvalue_expf[g,]))
#  nsurv[,g] = as.factor(ifelse(as.numeric(pvalue_expf[g,])>median_value,"High","Lower"))
#}

#nsurv2 = nsurv[,colnames(nsurv)!="patient"]
#multiCox=coxph(Surv(time, statu) ~ ., data = nsurv2)

#tdmultiCoxSum=summary(multiCox)
#outResult=data.frame()
#outResult=cbind(
#  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
#  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
#  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
#  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
#outResult=cbind(id=row.names(outResult),outResult)


pc <- prcomp(t(pvalue_expf))
value <- pc$x[, 1] * sign(cor(pc$x[, 1], colMeans(pvalue_expf,na.rm = T)))

result = as.data.frame(value)

result_final = data.frame(sample=rownames(result),result)
write.table(result_final,
            file="sample_PCA_score.txt",row.names = FALSE,sep="\t",
            quote = FALSE)


nsurv2 <- nsurv
#q1 <- quantile(nresult$value, 0.25)
#median_value <- quantile(nresult$value, 0.5)
#q3 <- quantile(nresult$value, 0.75)

nsurv2$group = ifelse(result$value>median(result$value),"High","Lower")

nsurv2$group = factor(nsurv2$group,levels = c("High","Lower"))

nsurv2$status <- grepl("1", nsurv2$statu, ignore.case = TRUE)
colnames(nsurv2) = c("statu","Patient","days_to_death","group","status")

fit <- survfit(Surv(days_to_death,status) ~ group,  data = nsurv2) 

pdf("survival_plot.pdf")

p = ggsurvplot(fit, # 创建的拟合对象
           title = "ImmunoScore Model",
           data = nsurv2,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           xlab="Time(days)",
           ylab = "Overall survival",
           pval = TRUE, # 添加P值
           risk.table = FALSE, # 绘制累计风险曲线
           surv.median.line = "hv", # 添加中位生存时间线
           add.all = FALSE, # 添加总患者生存曲线
           palette = "nejm",
           ggtheme=custom_theme(),
           legend = c(0.8,0.75),
           legend.labs=c("High","Low"),
           legend.title="Score")

print(p,newpage = FALSE)
dev.off()

#print(p)
#ggsave("survival_plot.pdf", width = 18, height = 18, units = "cm")

