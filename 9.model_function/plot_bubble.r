
rm(list=ls())
library(ggplot2)
library(stringr)

theme <- theme(
  panel.background=element_blank(),
  panel.border=element_rect(fill=NA),
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.x = element_text(colour="black",size=12,angle = 80,vjust = 0.4,hjust = 0.4),
  axis.text.y = element_text(colour="black",size=15),
  axis.ticks = element_line(colour="black"),
  # plot.margin = unit(c(1,1,1,1),"line"),
  # legend.title = element_blank(),
  # legend.key = element_rect(fill="white")
)

infile = "gProfiler_hsapiens_2024-1-31 16-32-33__intersections_selected.csv"
xname = "GeneRatio"

pathway = read.csv(infile,header=T,stringsAsFactors=FALSE)
pathway$GeneRatio = format(as.double(pathway$intersection_size/pathway$query_size),scientific = FALSE,digits = 2)

p2 = ggplot(pathway,aes(as.character(GeneRatio),reorder(term_name,as.numeric(GeneRatio))))
p=p2 + geom_point(aes(size=intersection_size,color=-1*log10(adjusted_p_value)))

pr = p+scale_color_gradient(low="blue",high = "red")
pr = pr+labs(color=expression(-log[10](adjusted_p_value)),
             size="Gene number",x=xname,y="")
pr=pr + theme_bw() + theme
png(width=2600,height=1300,"go_kegg_plot_bubble.png",res=300)
pr
dev.off()
tiff(width=2600,height=1300,"go_kegg_plot_bubble.tiff",res=300)
pr
dev.off()
pdf("go_kegg_plot_bubble.pdf",width = 10,height = 8)
pr
dev.off()
