##module load  R/3.6.2

library(Seurat)

library(Matrix)

library(tidyverse)

library(diem)

## Load the datasets
basefolder <- "../kallisto/bus/"
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


read_count_output <- function(dir, name) {
      dir <- normalizePath(dir, mustWork = TRUE)
        m <- readMM(paste0(dir, "/", name, ".mtx"))
        m <- Matrix::t(m)
        m <- as(m, "dgCMatrix")
        # The matrix read has cells in rows
        ge <- ".genes.txt"
        genes <- readLines(file(paste0(dir, "/", name, ge)))
        barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
        colnames(m) <- barcodes
        rownames(m) <- genes
        return(m)
      }


adata <- map(expNames,function(ii){
    ##
    expPrefix = ii;
    dir= paste0(folders[ii],"/counts_unfiltered")
    cat("#Loading ",dir,"... \n")
    spliced = read_count_output(dir,"spliced")
    ## unspliced = read_count_output(dir,"unspliced")
    sce <- create_SCE(spliced)
    cat(dim(sce),"\n")
    ## Remove debris...  New diem version seems different, we could try new or SoupX.
    sce <- diem(sce, top_n=16000)
    dim(sce)
    ##
    sc = convert_to_seurat(sce)
    cat("#Final: ",dim(sc),"\n")
    sc
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
