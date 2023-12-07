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

#load("./3_MergeDemux_Output/scFilteredSeurat.Rdata")

# soupx
sc<-read_rds("1_soupx/seuratObj-after-mt-filtering.2022-03-10.rds")
sc <- subset(sc1, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)
# DefaultAssay(sc1) <- "RNA"



outFolder="./4_harmony_res0.4/"
system(paste0("mkdir -p ", outFolder))

### Merge


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


#sapply(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),function(x)
#sapply(c(1,1.1,1.2),function(x)
#sapply(c(0.5,0.6,0.8,"1.0"),function(x)
sapply(c(1.5,2,2.5),function(x)
  {
  #outFolder=paste0("./4_harmony_res",x,"/")
  outFolder=paste0("./4_harmony_SoupX_res",x,"/")
  system(paste0("mkdir -p ", outFolder))
  
  sc2 <- FindClusters(sc, verbose = TRUE,resolution=as.numeric(x))
  fname=paste0(outFolder,"sc.NormByLocationRep.Harmony.rds")
  write_rds(sc2,fname)
  
  sc2 <- subset(sc2, subset = Labor%in% c("TIL" ,"TNL"))
  
  ## Make a simple plot here:
  fname=paste0(outFolder,"UMAP_LocationHarmony.pdf");
  pdf(fname,width=10,height=5)
  aa <- FetchData(sc2,c("UMAP_1","UMAP_2","seurat_clusters","Location","Labor","SNG.BEST.GUESS")) 
  
  p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
    geom_point(size=0.1) +
    facet_grid(Labor ~ Location) +
    theme_bw()
  p1
  ##    theme_black()
  dev.off()
  
  
})


################
sc <- FindClusters(sc, verbose = TRUE,resolution=0.4)
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