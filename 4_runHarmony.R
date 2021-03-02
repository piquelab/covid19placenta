######################################
### run Harmony - without filtering ###
######################################


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

load("../2020-10-02/3_MergeDemux_Output/scFilteredSeurat.Rdata")


outFolder="./4_harmony/"
system(paste0("mkdir -p ", outFolder))

### Merge

dim(sc)

table(sc$Library)

sc@meta.data$Location <- "CAM"
sc@meta.data$Location [grepl("PVBP",sc@meta.data$EXP)] <- "PVBP" 

table(sc$Library,sc$Location) 


## Harmony

DefaultAssay(sc) <- "RNA"

#LogNormalize
sc <- NormalizeData(sc, verbose=TRUE) 

#Identifies features that are outliers on a 'mean variability plot'.

#vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

#Scales and centers features in the dataset. If variables are provided in vars.to.regress, they are individually regressed against each feautre, and the resulting residuals are then scaled and centered.
sc <- ScaleData(sc, verbose = TRUE) 

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)


# remove the influence of dataset-of-origin from the embedding. By default, Harmony accepts a normalized gene expression matrix and performs PCA. Since here we already have the PCs, we specify do_pca=FALSE. The matrix harmony_embeddings is the matrix of Harmony corrected PCA embeddings.
sc <- RunHarmony(sc,c("Library"),reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)


###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE,resolution=0.6)

################

## save object
fname=paste0(outFolder,"sc.NormByLocationRep.Harmony.rds")
write_rds(sc,fname)


#################

## Make a simple plot here:
fname=paste0(outFolder,"UMAP_LocationHarmony.pdf");
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
md<-sc@meta.data
md<-md[,c("Library","nCount_RNA")]
md  %>% dplyr::count(Library,nCount_RNA )

md %>% 
  group_by(Library) %>% 
  summarise(UMI = sum(nCount_RNA))