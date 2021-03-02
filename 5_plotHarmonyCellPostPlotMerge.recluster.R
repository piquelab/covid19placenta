## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

######################################
### plot cell types ###
## using resolution 0.8
######################################


library("BiocParallel")
register(MulticoreParam(12))

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


sc <- read_rds("../2020-10-02/4_harmony/sc.NormByLibrary.Harmony,StringentFiltering.rds")


## Move this part to the runHarmony.R script
resolutions <- c(0.8,0.6,0.4,0.3,0.2,0.1)
## Default resolution is 0.8 K=28 , run also 0.6 K=24, 0.4 K=23, 0.2 K=14, 0.1 K=10, 0.3 K=20
sc <- FindClusters(sc, verbose = TRUE,resolution=0.8)


##pred.singler <- read_rds("./4_harmony_cellClass.v3/sc.NormByLibary.ref.Harmony_singler.rds") %>%
##    as.data.frame()

md <- read_rds("../2020-10-02/4_harmony_cellClass/sc.NormByLocation.ref.Harmony.singler.rds") %>%
    as.data.frame %>%
    rownames_to_column("BARCODES") %>%
    select(BARCODES,pruned.labels,scLabor_ID=pruned.labels,scLabor_Score=tuning.scores.second)

md <- sc@meta.data %>% rownames_to_column("BARCODES") %>%
    left_join(md) 

identical(md$BARCODES,rownames(sc@meta.data))


outFolder="./5_harmony_cellClass_plots_res0.8/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

dim(sc)

table(sc$Library)

table(sc$Location) 

table(md$scLabor_ID, md$Location)

table(md$scLabor_ID, md$Origin)

table(md$scLabor_ID, md$Location)

table(md$SNG.BEST.GUESS, md$scLabor_ID)




cc <- md %>% select(seurat_clusters,scLabor_ID) %>%
    group_by(seurat_clusters,scLabor_ID) %>%
    summarize(n=n()) %>%
    group_by(seurat_clusters) %>%
    mutate(prop=n/sum(n)) %>%
    arrange(-n)

# majority of ce
cc2 <- cc %>% top_n(n=1) %>% mutate(cname=paste0(seurat_clusters,":",scLabor_ID)) 
clust2name <- cc2$cname
names(clust2name)=cc2$seurat_clusters
clust2name

source("./scLaborColors.R")

cc2 %>% select(seurat_clusters,scLabor_ID) %>%
    arrange(seurat_clusters) %>%
    mutate(color=group.colors[scLabor_ID]) %>%
    write_tsv(paste0(outFolder,"ClusterAssignment.tsv"))

cc2<-cc2 %>% select(seurat_clusters,scLabor_ID) %>%
    arrange(seurat_clusters) %>%
    mutate(color=group.colors[scLabor_ID])


sc@meta.data <- md
write_rds(sc,paste0(outFolder,"SeuratObject.rds"))

tt <- table(md$seurat_clusters,md$scLabor_ID)
tt <- tt/rowSums(tt)

library(pheatmap)

fname=paste0(outFolder,"Singler.HeatMap.pdf");
pdf(fname,width=7,height=4)
pheatmap(t(tt),cluster_rows=TRUE,cluster_cols=FALSE,scale="none")
dev.off()

##pbmc <- RenameIdents(pbmc, new.cluster.ids)
fname=paste0(outFolder,"UMAP_Harmony.png");
png(fname,width=1000,height=1000)
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()



## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.png");
png(fname,width=1600,height=1200)
aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","Origin","status")) 
aa$seurat_clusters <- clust2name[aa$seurat_clusters]
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Location)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Location")) +
##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_ConditionHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Condition)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()

## Make a simple plot here:
fname=paste0(outFolder,"UMAP_OriginHarmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Origin)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Origin")) +
##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()


## Make a simple plot here:
fname=paste0(outFolder,"UMAP_StatusSoC_Harmony.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=status)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="status")) +
##    facet_wrap(~LocTime) +
    theme_bw()
p1
##    theme_black()
dev.off()



## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.Cell_annotation.png");
png(fname,width=1600,height=1200)
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
    geom_point(size=0.1) +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell type")) +
    facet_wrap(Condition ~ Location) +
    theme_bw()
p1
##    theme_black()
dev.off()



fname=paste0(outFolder,"UMAP_Location.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=seurat_clusters,fill=Origin)) +
    geom_bar(position="stack") +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell origin")) +
    facet_grid(.~Location) + coord_flip() +
    theme_bw()
p2
##    theme_black()
dev.off()


fname=paste0(outFolder,"UMAP_LocationCondition.Barplot.pdf");
pdf(fname,width=10,height=6)
p2 <- ggplot(aa,aes(x=seurat_clusters,fill=Condition)) +
    geom_bar(position="stack") +
##    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Condition")) +
    facet_grid(.~Location) + coord_flip() +
    theme_bw()
p2
##    theme_black()
dev.off()





### END- HERE ###
########################################################


