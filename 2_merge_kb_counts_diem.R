library(Seurat)

library(Matrix)

library(tidyverse)

library(diem)

## Load the datasets
basefolder <- "/nfs/rprdata/scilab/labor2/kallisto.covid19/bus/"
expNames <- dir(basefolder,"^C19C*")


expNames

subfolder <- ""
folders <- paste0(basefolder,expNames,subfolder)
ind <- file.info(folders)$isdir
ind[is.na(ind)]<- FALSE
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames
folders

adata <- map(expNames,function(ii){
    ##
    expPrefix = ii;
    folders[expNames[1]]
    fName= paste0(folders[ii],"/counts_unfiltered/adata.h5ad")
    cat("#Loading ",fName,"... \n")
    adata = ReadH5AD(file = fName,verbose=TRUE)
    cat(dim(adata),"\n")
    sce <- create_SCE(rbind(adata@assays$spliced@data,adata@assays$unspliced@data))
    cat(dim(sce),"\n")
    ## Remove debris...  And this seems to have changed ... or consider switching to SoupX. 
    sce <- diem(sce)
    ##
    sc = convert_to_seurat(sce)
    cat("#Final: ",dim(sc),"\n")
    adata <- adata[,colnames(sc)]
    adata
})


names(adata) <- expNames


system("mkdir -p 2_kb_diem_Output/")
write_rds(adata,"2_kb_diem_Output/kb_diem_Seurat.list.rds")

sparse.size <- object.size(adata)
sparse.size

sc <- merge(adata[[1]], y = adata[-1], project = "KallistoDiem",add.cell.ids=expNames)

dim(sc)

system("mkdir -p 2_kb_diem_Output/")
write_rds(sc,"2_kb_diem_Output/kb_diem_Seurat.obj.rds")

##############################
