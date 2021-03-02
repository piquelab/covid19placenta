###################################################
### cell type plots- finalizing the cell lables ###
###################################################

## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)

#################
##library(SingleR)



future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################
## Testing sc transform           
## 2and3_Diem_Output
## adata <- read_rds("./kb_diem_Output/kb_diem_Seurat.list.rds")

sc <- read_rds("../2020-10-02/5_harmony_cellClass_plots_res0.8/SeuratObject.rds")
identical(colnames(sc),rownames(sc@meta.data))


new_names <- read_tsv("../2020-10-02/5_harmony_cellClass_plots_res0.8/ClusterAssignment.res0.8.tsv")
clust2Names <- new_names$scLabor_ID
names(clust2Names) <- new_names$seurat_clusters
cc <- new_names %>% select(scLabor_ID,color) %>% unique 
cluster.Colors <- cc$color
names(cluster.Colors) <- cc$scLabor_ID

tempcol<-cluster.Colors["B-cell"]
cluster.Colors["B-cell"]<-cluster.Colors["T-cell"]
cluster.Colors["T-cell"]<-tempcol

sc@meta.data$cluster_name <- clust2Names[sc@meta.data$seurat_clusters]

table(sc$Library)

table(sc$Location) 

table(sc$cluster_name, sc$Location)

table(sc$cluster_name, sc$Origin)

table(sc$SNG.BEST.GUESS, sc$cluster_name)


outFolder="./6_harmony_rename_res0.8_plots/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

write_rds(sc,paste0(outFolder,"SeuratObject.rds"))

Idents(sc) <- "cluster_name"

fname=paste0(outFolder,"UMAP_Harmony.png");
png(fname,width=1000,height=1000)
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

aa <- FetchData(sc,c("UMAP_1","UMAP_2","Location","Condition","Origin","status","FetalSex","seurat_cluster","cluster_name")) 
head(aa)


fname=paste0(outFolder,"UMAP_LocationHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
    geom_point(size=0.1) +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    theme_bw()
p1
##    theme_black()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
    geom_point(size=0.1) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    theme_bw()
p1
##    theme_black()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
    geom_point(size=0.1) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    scale_color_manual(values=c("Control"="#333399","w/COVID-19"="#A50021"))

p1
##    theme_black()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_OriginHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Origin)) +
    geom_point(size=0.1) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Origin")) +
    theme_bw()
p1
##    theme_black()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_StatusSoC_Harmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=status)) +
    geom_point(size=0.1) +
    guides(colour = guide_legend(override.aes = list(size=5),title="status")) +
    theme_bw()
p1
##    theme_black()
dev.off()



## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.Cell_annotation.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=cluster_name)) +
    geom_point(size=0.1) +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell type")) +
    facet_grid(Condition ~ Location) +
    theme_bw()
p1
##    theme_black()
dev.off()



fname=paste0(outFolder,"UMAP_Location.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=cluster_name,fill=Origin)) +
    geom_bar(position="stack") +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell origin")) +
    facet_grid(.~Location) + coord_flip() +
    theme_bw()
p2
##    theme_black()
dev.off()


 
fname=paste0(outFolder,"UMAP_LocationCondition.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=cluster_name,fill=Condition)) +
    geom_bar(position = position_stack(reverse = TRUE)) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    facet_grid(.~Location) + coord_flip() +
    scale_fill_manual(values=c("Control"="#333399","w/COVID-19"="#A50021"))
p2
##    theme_black()
dev.off()

### END- HERE ###
########################################################


