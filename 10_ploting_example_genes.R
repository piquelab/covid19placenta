##################################################################
### Ploting some top 3 down-regulated and up-regulated genes ###
##  Plotsing DE genes per cell types
##################################################################

library(Seurat)
library(Matrix)
library(tidyverse)

library(tidyverse)
library(dplyr)
library(stringr)

new_names <- read_tsv("../2020-10-02/5_harmony_cellClass_plots_res0.8/ClusterAssignment.res0.8.tsv")
clust2Names <- new_names$scLabor_ID
names(clust2Names) <- new_names$seurat_clusters
cc <- new_names %>% select(scLabor_ID,color) %>% unique 
cluster.Colors <- cc$color
names(cluster.Colors) <- cc$scLabor_ID
names(cluster.Colors) <- cc$scLabor_ID
tempcol<-cluster.Colors["B-cell"]
cluster.Colors["B-cell"]<-cluster.Colors["T-cell"]
cluster.Colors["T-cell"]<-tempcol

location.colors  = c( "PVBP"="#fc8d62", "CAM"="#8da0cb")
condition.colors = c("Control"="#333399","w/COVID-19"="#4daf4a")




filter_sample<-"HPL20874"
outFolder <- paste0("./8_example_genes_Plots_",paste(filter_sample,collapse="_"),"/")

system(paste0("mkdir -p ",outFolder))


##loading all DE genes
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_filtered_HPL20874_2020-11-28/SIG.combined.2020-11-28.tsv")


# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
res<-res[order(abs(res$log2FoldChange), abs(res$baseMean),decreasing = TRUE),]



# choosing 3 up-regulated and 3 down-regulated genes from res (all combined DE genes from categories)
selected_cell_types<-c("T-cell","Macrophage-2","Macrophage-1")
res_examplegenes<-res  %>% filter (Cell_type %in% selected_cell_types)


t_cell_genes<-c("THBS1","LRATD2","CMSS1","FARSA","PTRH1","PMS2P1")
m2_genes<-c("C1orf21","DUSP4","SMOX","BHLHE41","TXNIP","TNFAIP8L2")
m1_genes<-c("CMSS1","NCS1","TRAF5","MNDA","BHLHE41" ,"CLEC4E" )


res_t_cell<-res %>% filter (Cell_type=="T-cell"  & gene_name %in% t_cell_genes )
res_t_cell<-res_t_cell[order(res_t_cell$log2FoldChange),]
res_m1<-res %>% filter (Cell_type=="Macrophage-1"  & gene_name %in% m1_genes )
res_m1<-res_m1[order(res_m1$log2FoldChange),]
res_m2<-res %>% filter (Cell_type=="Macrophage-2"  & gene_name %in% m2_genes )
res_m2<-res_m2[order(res_m2$log2FoldChange),]

res_violin<-rbind(res_t_cell,res_m2,res_m1)

## Load cells
sc <- read_rds("../2020-10-02/6_harmony_rename_res0.8_plots/SeuratObject.rds")

#meta data
md <- sc@meta.data




##
a2 <- FetchData(sc,c("UMAP_1","UMAP_2","cluster_name","Location","Condition","Origin","Pregnancy_ID","Condition")) 

myscale = 1/colSums(sc@assays$RNA@counts)*1000000



##rec <- map_dfr(1:nrow(res),function(ii){
rec <- map_dfr(1:(dim(res_violin)[1]),function(ii){
 a2 <- FetchData(sc,c("UMAP_1","UMAP_2","cluster_name","Location","Condition","Origin","Pregnancy_ID")) 
 a2$gene_name=res_violin$gene_name[ii]
 a2$cname=paste(res_violin$gene_name[ii],res_violin$cname[ii],sep="_")
 a2$Expression=sc@assays$RNA@counts[res_violin$kbid[ii],]*myscale
 a2$log2FoldChange<-res_violin$log2FoldChange[ii]
 a2$baseMean<-res_violin$baseMean[ii]
 aux <- a2 %>% filter(Location==res_violin$Location[ii],Origin==res_violin$Origin[ii],cluster_name==res_violin$Cell_type[ii])
## a2$Expression=sc2@assays$RNA@scale.data[selGenes$kbid[ii],]
 aux
})
dim(rec)



## to make the boxplot


######################################################
### boxplot
######################################################



pdf(paste0(outFolder,"boxplot_Top10DE.pdf"),width=25,height=10)
#gg <-
rec %>% ggplot(aes(x=cname,y=log10(Expression),fill=Condition)) + 
geom_boxplot()+
scale_fill_manual(values=c("Control"="#333399","w/COVID-19"="#A50021"))+
theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
theme(axis.title.x=element_blank())
dev.off()

######################################################
### splited violin pot
######################################################

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                          draw_group = function(self, data, ..., draw_quantiles = NULL) {
                            data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                            grp <- data[1, "group"]
                            newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                            newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                            newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                            
                            if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                              stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                        1))
                              quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                              aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                              aesthetics$alpha <- rep(1, nrow(quantiles))
                              both <- cbind(quantiles, aesthetics)
                              quantile_grob <- GeomPath$draw_panel(both, ...)
                              ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                            }
                            else {
                              ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                            }
                          })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



#rec<-rec[order(rec$cluster_name , decreasing = TRUE),]
#rec<-rec  %>% filter (cluster_name %in% c("T-cell","Macrophage-1","Macrophage-2"))

rec<-rec %>% group_by(desc(cluster_name)) %>% arrange(log2FoldChange,baseMean, .by_group = TRUE) 


rec$cluster_name <- factor(rec$cluster_name,levels=unique(rec$cluster_name))

pdf(paste0(outFolder,"violinplot_genes_per_celltypes_v2.pdf"),width=25,height=10)
rec %>% ggplot(aes(x=reorder(gene_name,(log2FoldChange)),y=log10(Expression),fill=Condition)) + 
  geom_split_violin()+
  scale_fill_manual(values=c("Control"="#333399","w/COVID-19"="#A50021"))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +  #,text = element_text(size=30)
  facet_grid(~cluster_name, scale="free",space="free") 
dev.off()



unique(rec[(rec$cluster_name=="Macrophage-1"),c("gene_name","log2FoldChange")])
unique(rec[(rec$cluster_name=="Macrophage-2"),c("gene_name","log2FoldChange")])
unique(rec[(rec$cluster_name=="T-cell"),c("gene_name","log2FoldChange")])


gg <- rec %>% ggplot(aes(x=reorder(Pregnancy_ID,-log10(Expression)),y=log10(Expression),fill=Condition)) +
    geom_boxplot(position=position_dodge(1)) +
    scale_fill_manual(values=c("Control"="#333399","w/COVID-19"="#A50021"))+
    theme(axis.text.x = element_text(angle = 90))+
    facet_wrap( ~ cname,scales = "free",ncol=2)+
    #scale_x_discrete(expand=c(0.0001,0.01))
    scale_x_discrete(expand=c(0, 0.5))+
    geom_abline(0,1) +
    theme(text = element_text(size = 7))
theme(panel.spacing =unit(.1, "lines"))
#theme(panel.spacing.x=unit(1, "lines"),panel.spacing.y=unit(1, "lines"))

ggsave(paste0(outFolder,"Top10DEplot_WithPregnancyID.png"),gg)



### THIS is to make a summary per gene_location. 
sum_rec <- rec %>% group_by(cname,Condition) %>% summarize(Prop=mean(Expression>0),Expr=mean(Expression[Expression>0]))

dim(sum_rec)

head(sum_rec)

#sum_rec <- sum_rec %>% left_join(selGenes) 
write_csv(sum_rec,paste0(outFolder,"top10DE_summary.csv"))


######################################################
## forestPlot
######################################################
DE_number_per_celltype<-tapply(res$padj ,factor(res$Cell_type),length )
DE_number_per_celltype["Monocyte"]<-12
res$DE_number_per_celltype<-DE_number_per_celltype[res$Cell_type]

res<-res[order(res$DE_number_per_celltype,decreasing = TRUE),]


res$Cell_type <- factor(res$Cell_type,levels=unique(res$Cell_type))

#coloring
p<-res %>% 
  ggplot(aes(x=reorder(gene_name,(log2FoldChange)),y=log2FoldChange,
             alpha=(padj<0.1))) +
  geom_point(aes(size=padj),position=position_dodge(width=0.2)) +
  scale_size("q-values", trans="log10", range=c(4, 0.2),limits=c(1E-10,1), breaks=c(1E-12,1E-6,0.001,0.01,0.1)) +
  scale_color_manual(guide = guide_legend(reverse = TRUE) ) +
  geom_errorbar(aes(ymax = log2FoldChange + 1.96*lfcSE, ymin = log2FoldChange - 1.96*lfcSE),width=0,position=position_dodge(width=0.2)) +
  scale_alpha_manual(values=c(0.3, 1.0),guide=FALSE) +
  geom_hline(yintercept=0,lty=2) + 
  xlab("Gene") + ylab(expression(log[2](Fold~Change))) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0,hjust=0,vjust=0.5),strip.background=element_rect(fill="white",color="white")) +
  coord_flip() + 
  facet_grid(Cell_type ~ . ,scales = "free_y",space="free") 
ggsave(paste0(outFolder,"forestPlot_DEgenes_cname_v2.pdf"),p,width=7,height=14)

g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("T-cell"="red","Macrophage-2"="tomato","Macrophage-1"="hotpink","Monocyt"="lightslateblue","LED"="yellow4","CTB"="purple4" , "Stromal-3"="peru","npiCTB"="limegreen" )
k <- 1

for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col<-fills[k]
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$lwd=4
  k <- k+1
}
ggsave(paste0(outFolder,"forestPlot_DEgenes_cname_v2.pdf"),g,width=7,height=14)

