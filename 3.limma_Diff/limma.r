#!/usr/bin/Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(Biobase))
suppressMessages(library(genefilter))
suppressMessages(library(limma))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
colors = colorRampPalette(brewer.pal(9,'GnBu'))(100)


up_color = "red"
down_color = "green"
not_color = "black"
## load args

setwd("E:/workspace/keyan/GC_tme/Article_results/3.Limma_screen_biomarker/")
exp_file="TCGA-STAD.htseq_fpkm.tsv_format.txt_tumor_mRNA_exp.txt"
pheno_file="immune_cluster_ordered.txt"
outdir="."
fc=1.5
pval=0.05
fc = log2(as.numeric(fc))
tcontrasts=c("Cluster1~Cluster2","Cluster1~Cluster3",
            "Cluster1~Cluster4","Cluster2~Cluster3",
            "Cluster2~Cluster4","Cluster3~Cluster4")
logmethod = "1"
prefix="immune_cluster"
xlim="0"
ylim="0"

## read and check expdata and pheodata
expr_data =  try(read.table(exp_file, header=T, quote="",
                             row.names=1, sep='\t',check.names=F,
                             as.is=TRUE))

pheno_data = try(read.table(pheno_file,check.names=F,
                            header=T,sep='\t',row.names=1,
                            as.is=TRUE))
colnames(pheno_data) = "group"

for (tcontrast in tcontrasts){
        noutdir = paste(outdir,tcontrast,sep="/")
        createdir = sprintf("mkdir -p %s",noutdir)
        system(createdir)
        con = unlist(strsplit(tcontrast,'~'))
        print(con)
        stats = c("Total","Diff","up","down")
        headers = pheno_data[,1] == con[1] | pheno_data[,1]==con[2]
        npheno_data = subset(pheno_data,group==con[1]|group==con[2])
        # phenodata = pheno_data[pheno_data$group == "control" | pheno_data$group == "model",]
        nexpr_data = expr_data[,as.character(rownames(npheno_data))]
        nexpr_data = as.matrix(nexpr_data)
        write.table(cbind(rownames(nexpr_data),nexpr_data),paste(noutdir,paste(prefix,"_rawdata.xls",sep=""),sep="/"),
            sep="\t",row.names=FALSE, quote=F)
        qx <- as.numeric(quantile(nexpr_data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) ||
                  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
                  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
        if (LogC & logmethod=="auto"){ 
          nexpr_data[which(nexpr_data <= 0)] <- NaN
          nexpr_data <- log2(nexpr_data)
          # expr_data[which(expr_data<=0)] <- NaN
          # expr_data <- log2(expr_data)
          print("auto transform log")
           }
        if (logmethod == "0"){
          print("force transform log")
            minv = min(expr_data)
            if (minv == 0){
                nexpr_data = log2(nexpr_data+0.001)
            }else{
                nexpr_data = log2(nexpr_data) 
            }
        }
        nexpr_data = na.omit(nexpr_data)
        # expr_data = na.omit(expr_data)
        normf = paste(noutdir,paste(prefix,"_Normalization.xls",sep=""),sep="/")
        write.table(cbind(rownames(nexpr_data),nexpr_data),normf,
            sep="\t",row.names=FALSE, quote=F)
        # write.table(cbind(id=rownames(expr_data),expr_data),paste(noutdir,paste(prefix,"_total_Normalization.xls",sep=""),sep="/"),
        #     sep="\t",row.names=FALSE, quote=F)
        lev = unique(npheno_data$group)
        g = factor(npheno_data$group,levels=con)
        design <- model.matrix(~0+g)
        colnames(design) = con
        fit <- lmFit(nexpr_data, design)
        contrast = paste(con,collapse = '-')
        contrast.matrix <- makeContrasts(contrast,
                                             levels=design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        tt = topTable(fit2, coef=1, n=Inf,
                          sort.by = 'logFC',adjust="BH")
        print(colnames(tt))
        tt <- subset(tt, select=c("logFC","P.Value","adj.P.Val"))
        diffa <- na.omit(tt)
        diff = tt[abs(tt$logFC) > fc & tt$adj.P.Val < pval,]
        diff = na.omit(diff)
        diff2 = tt[tt$logFC > fc & tt$adj.P.Val < pval,]
        diff2 = na.omit(diff2)
        diff3 = tt[tt$logFC < -fc & tt$adj.P.Val < pval,]
        diff3 = na.omit(diff3)
        outpatha = paste(noutdir,
                        paste(prefix,'_',tcontrast,'_all.xls',sep=''),
                              sep='/')
        outpath = paste(noutdir,
                        paste(prefix,'_',tcontrast,'_diff.xls',sep=''),
                              sep='/')
        outpath2 = paste(noutdir,
        				paste(prefix,"_",tcontrast,"_up.xls",sep=''),
        					  sep="/")
        outpath3 = paste(noutdir,paste(prefix,"_",tcontrast,"_down.xls",sep=''),sep="/")
        # print(outpath)
        # print(outpath2)
        # print(outpath3)
        statsf = paste(noutdir,paste(prefix,"_",tcontrast,"_stats.xls",sep=''),sep="/")
        diffaout = cbind(rownames(diff),diff)

        diffouta = cbind(genesymbol=rownames(diffa),diffa,fold=2**diffa$logFC)
        colnames(diffouta)[1] = "genesymbol" 


        diffout = cbind(genesymbol=rownames(diff),diff,fold=2**diff$logFC)
        colnames(diffout)[1] = "genesymbol" 

        diff2out = cbind(genesymbol=rownames(diff2),diff2,fold=2**diff2$logFC)
        colnames(diff2out)[1] = "genesymbol"

        diff3out = cbind(genesymbol=rownames(diff3),diff3,fold=2**diff3$logFC)
        colnames(diff3out)[1] = "genesymbol"

        stats = rbind(stats,c(nrow(expr_data),nrow(diff),nrow(diff2),nrow(diff3)))

        write.table(diff2out,outpath2,sep='\t',row.names=FALSE, quote=F)
        write.table(diff3out,outpath3,sep='\t',row.names=FALSE, quote=F)
        # diff$genenames = rownames(diff)
        write.table(diffout,outpath,sep='\t',row.names=FALSE, quote=F)
        write.table(diffouta,outpatha,sep="\t",row.names=FALSE,quote=F)


        write.table(stats,statsf,sep="\t",row.names=FALSE,quote=F,col.names=FALSE)
        matrix_file = paste(noutdir,
                          paste(prefix,'_',tcontrast,'_matrix.xls',sep=''),
                          sep='/')
        matrix_data = nexpr_data[rownames(diff),]
        matrix_data_out = cbind(rownames(matrix_data),matrix_data)
        colnames(matrix_data_out)[1] = "genesymbol"
        write.table(matrix_data_out,matrix_file,sep='\t',row.names=FALSE, quote=F)
        new_pheno_file = paste(noutdir,
                              paste(prefix,'_',tcontrast,'_pheno.xls',sep=''),
                              sep='/')
        new_pheno = pheno_data[headers,,drop=F]
            write.table(new_pheno,new_pheno_file,sep='\t',col.names=NA, quote=F)

            # mask <- with(tt, abs(logFC) > fc & P.Value < pval)
            # cols <- ifelse(mask,"red","black")
            cols <- ifelse((abs(tt$logFC)>fc & tt$adj.P.Val<pval),
                    ifelse(tt$logFC>0,up_color,down_color),not_color)
            png(paste(noutdir,
                      paste(prefix,'_',tcontrast,'_vocalno.png',sep=''),
                      sep='/'),width=480,height=640)

                max_p = max(-log10(na.omit(tt$adj.P.Val))) + 1 
                max_d = max(abs(na.omit(tt$logFC))) + 0.5
                if(max_p>12){
                    max_p = 12
                }
                  if(xlim != "0"){
                    max_d=as.numeric(xlim)
                }
                   if(ylim != "0"){
                    max_p=as.numeric(ylim)
                }
            par(mar=c(5, 5, 4, 2) + 0.1)
            p = plot(tt$logFC, -log10(tt$adj.P.Val), cex=.8, pch=16,
                          xlim=c(-max_d,max_d), ylim=c(0,max_p),
                          xlab="log2 Fold Change",
                          ylab="-log10 P-value",
                          main=tcontrast,ann = F,
                          col=cols,xaxt = "n", yaxt = "n")
            p = p + title(xlab= "log2 Fold Change", ylab ="-log10 P-value", main=tcontrast,cex.lab = 1.8,cex.main = 1.8,
                font.main= 1.8, col.main= "black")
            p = p + axis(1,cex.axis = 1.8)
            p = p + axis(2,cex.axis = 1.8)
            p = p + abline(h=-log10(pval),v=c(-fc,fc), lty=2)
            legend("topright",c("Up","Not","Down"),pch=16,col=c(up_color,not_color,down_color),bty="n",cex=1.8,xpd=T)
            print(p)
            dev.off()
            
            pdf(paste(noutdir,
                       paste(prefix,'_',tcontrast,'_vocalno.pdf',sep=''),
                       sep='/'),width=7.5,height=10)
            par(mar=c(5, 5, 4, 2) + 0.1)
            p = plot(tt$logFC, -log10(tt$adj.P.Val), cex=.8, pch=16,
                          xlim=c(-max_d,max_d), ylim=c(0,max_p),
                          xlab="log2 Fold Change",
                          ylab="-log10 P-value",
                          main=tcontrast,ann = F,
                          col=cols,xaxt = "n", yaxt = "n")
            p = p + title(xlab= "log2 Fold Change", ylab ="-log10 P-value", main=tcontrast,cex.lab = 1.8,cex.main = 1.8,
                font.main= 1.8, col.main= "black")
            p = p + axis(1,cex.axis = 1.8)
            p = p + axis(2,cex.axis = 1.8)
            p = p + abline(h=-log10(pval),v=c(-fc,fc), lty=2)
            legend("topright",c("Up","Not","Down"),pch=16,col=c(up_color,not_color,down_color),bty="n",cex=1.8,xpd=T)
            print(p)
            dev.off()


}






