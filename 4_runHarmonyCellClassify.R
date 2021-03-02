##############################################
### cell type classification using SingleR ###
### reference: parturition elife           ### 
##############################################

## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

library(harmony)

#################
library(SingleR)



future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################
## Testing sc transform           
## 2and3_Diem_Output
## adata <- read_rds("./kb_diem_Output/kb_diem_Seurat.list.rds")

load("../2020-10-02/3_MergeDemux_Output/scFilteredSeurat.Rdata")



sc3 <- read_rds("/nfs/rprdata/scilab/novogene/Analyses/Roger_20200218/3_scTransferLabel_scLabor/ST_Integrated.scLabor.obj.rds")

# Subset a Seurat object
sc3 <- subset(sc3, subset = nFeature_RNA > 100)



outFolder="./4_harmony_cellClass/"
system(paste0("mkdir -p ", outFolder))

### Merge
# to identify the cell types of sc1, another study with known celltypes will be merged to this study and the cell types will be identified

sc<-read_rds("../2020-10-02/4_harmony/sc.NormByLocationRep.Harmony.rds")
sc1 <- sc

##sc <- merge(sc1,list(sc2,sc3,sc4))

sc <- merge(sc1,list(sc3))


dim(sc)

table(sc$Library)

## table(sc$Location) 

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

he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

#unknown cell types
query.he <- he[,is.na(sc@meta.data$FinalName)]

#known cell types
ref.he <- he[,!is.na(sc@meta.data$FinalName)]

ref.labels <- sc@meta.data$FinalName[!is.na(sc@meta.data$FinalName)]

pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.rds")
write_rds(pred.labels,fname)

md <- pred.labels %>% as.data.frame() %>% 
    rownames_to_column("BARCODES") %>%
    left_join(sc@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"sc.NormByLocation.ref.Harmony.singler.csv")
write_csv(md,fname)

## save object.

sc<-read_rds("../2020-10-02/4_harmony_cellClass/sc.NormIntegrated.ref.Harmony.rds")

fname=paste0(outFolder,"sc.NormIntegrated.ref.Harmony.rds")
write_rds(sc,fname)


