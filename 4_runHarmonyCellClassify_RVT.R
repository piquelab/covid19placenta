#############################################################
### cell type classification using SingleR 
### reference: Vento-Tormo et al. Single-cell reconstruction of the early maternal-fetal interface in humans. Nature 563, 347-353 (2018).
#############################################################
library(Seurat)
library(Matrix)
library(tidyverse)
library(future)
library(harmony)
library(SingleR)

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

#RVT: Vento-Tormo, Roser, Mirjana Efremova, Rachel A. Botting, Margherita Y. Turco, Miquel Vento-Tormo, Kerstin B. Meyer, Jong-Eun Park et al. "Single-cell reconstruction of the early maternal–fetal interface in humans." Nature 563, no. 7731 (2018): 347-353

###########################################
#outFolder="./4_harmony_cellClass_RVT/"
#system(paste0("mkdir -p ", outFolder))

          
# reference
sc2 <- read_rds("/nfs/rprdata/scilab/novogene/otherdata/roser/Analysis/20200226/3_scTransferLabel_RoserPrep/ST_Integrated_RVT.obj.rds")
sc2 <- subset(sc2, subset = nFeature_RNA > 200)


# query: sc 
# load("3_MergeDemux_Output/scFilteredSeurat.Rdata")
# sc1 <- sc


# soupx
#sc1<-read_rds("1_soupx/seuratObj-after-mt-filtering.2022-03-10.rds")
#sc1<-read_rds("1_soupx/seuratObj-newfilter-merge.2022-03-10.rds")
#sc1 <- subset(sc1, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 & DIFF.LLK.BEST.NEXT > 3 & percent.mt < 25)

load("/wsu/home/groups/prbgenomics/covid19/covid19analysis_public_repo/3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc_covid<-sc
sc_covid <- subset(sc_covid, subset = Condition=="Control")
sc_covid@meta.data$Location <- "CAM"
sc_covid@meta.data$Location [grepl("PVBP",sc_covid@meta.data$EXP)] <- "PVBP" 
load("./3_MergeDemux_Output/scFilteredSeurat.Rdata")
sc@meta.data$Condition<-sc$Labor
sc1 <- merge(sc_covid,sc, project="parturition")


### Merge
sc <- merge(sc1,list(sc2))
dim(sc)
table(sc$Library)


## Harmony

DefaultAssay(sc) <- "RNA"

sc <- NormalizeData(sc, verbose=TRUE) 

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE) 

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

##sc <- RunHarmony(sc,c("Location","percent.mt","Rep"),reduction="pca")

sc <- RunHarmony(sc,c("Library"),reduction="pca",seed=TRUE)

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)


###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)




#sapply(c(0.5,0.6,0.8, 1.0,1.1),function(x)
#sapply(c(0.6,0.8,1.0,1.5),function(x)
sapply(c( 0.8, "1.0", 1.5,2),function(x)
{
  
#outFolder=paste0("./4_harmony_cellClass_RVT",x,"/")

#outFolder=paste0("./4_harmony_cellClass_SoupX_RVT",x,"/")


outFolder=paste0("./4_harmony_cellClass__with_covidcontrol_RVT",x,"/")

system(paste0("mkdir -p ", outFolder))
  

sc <- FindClusters(sc, verbose = TRUE, resolution=as.numeric(x))

################

he <- t(sc@reductions$harmony@cell.embeddings[,1:30])

query.he <- he[,is.na(sc@meta.data$annotation)]

ref.he <- he[,!is.na(sc@meta.data$annotation)]

ref.labels <- sc@meta.data$annotation[!is.na(sc@meta.data$annotation)]

pred.labels <- SingleR(test = query.he, ref = ref.he, labels = ref.labels)

##table(pred.labels)

table(pred.labels$pruned.labels)

sum(is.na(pred.labels$pruned.labels))

fname=paste0(outFolder,"sc.NormByLibrary.ref.Harmony.singler.RVT.rds")
write_rds(pred.labels,fname)

md <- pred.labels %>% as.data.frame() %>% 
    rownames_to_column("BARCODES") %>%
    left_join(sc1@meta.data %>% rownames_to_column("BARCODES"))


fname=paste0(outFolder,"sc.NormByLibrary.refRVT.Harmony.singler.csv")
write_csv(md,fname)

## save object.
fname=paste0(outFolder,"sc.NormByLibFullIntegrated.refRVT.Harmony.rds")
write_rds(sc,fname)
})



