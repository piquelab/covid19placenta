
######################################
### run Harmony - filtering doublet ###
######################################

## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

###########################################
## Testing sc transform           
## 2and3_Diem_Output
## adata <- read_rds("./kb_diem_Output/kb_diem_Seurat.list.rds")

load("3_MergeDemux_Output/scFilteredSeurat.Rdata")

outFolder="./4_harmony/"
system(paste0("mkdir -p ", outFolder))
##setwd(outFolder)

### Merge

sc <- subset(sc, subset = nFeature_RNA > 400 & percent.mt < 10 & status=="singlet")

dim(sc)

table(sc$Library)

sc@meta.data$Location <- "CAM"
sc@meta.data$Location [grepl("PVBP",sc@meta.data$EXP)] <- "PVBP" 

table(sc$Library,sc$Location) 

##table(sc$sclabor.tlabel)

##table(sc$Location,sc$sclabor.tlabel)

## Harmony

DefaultAssay(sc) <- "RNA"

sc <- NormalizeData(sc, verbose=TRUE) 

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE) 

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

##sc <- RunHarmony(sc,c("Location","percent.mt","Rep"),reduction="pca")
sc <- RunHarmony(sc,c("Library"),reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)


###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE,resolution=0.6)

################

## save object
fname=paste0(outFolder,"sc.NormByLibrary.Harmony,StringentFiltering.rds")
write_rds(sc,fname)


#################

## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LibraryHarmonyMoreStringent.pdf");
pdf(fname,width=12,height=5)
aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location","Condition","SNG.BEST.GUESS")) 
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
    geom_point(size=0.1) +
    facet_grid(Condition ~ Location) +
    theme_bw()
p1
##    theme_black()
dev.off()


### END- HERE ###
########################################################


