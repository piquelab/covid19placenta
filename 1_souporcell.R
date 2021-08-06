######################################
### souporcell for downstream analysis ###
######################################

##library("rhdf5")
library(Matrix)
library(tidyverse)
library(parallel)
library(data.table)

rm(list=ls())

outdir <- "./1_souporcell_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


#################################
### 1. folders of counts data ###
#################################
basefolder <- "../counts_cellranger/souporcell/"
demuxfn <- list.files(basefolder,pattern="^C19C*",include.dirs=TRUE)

expNames<-demuxfn
 
############################
### 2, read souporcell data ###
############################
demux <- mclapply(expNames,function(ii){
    cat("#Loading ", ii, "\n")
    fn <- paste0(basefolder, ii,"/clusters.tsv")
    dd <- data.frame(fread(fn,header=T))
    dd <- dd %>%
      mutate(barcode=paste0(ii,"_", gsub("-1","",barcode)),EXP=ii) %>%
      mutate(assignment=paste0(ii,"_",assignment)) %>%
      select(barcode,status,assignment,log_prob_singleton,log_prob_doublet)             
    ##
    fn <- paste0("cat ",basefolder, ii,"/plink2.kin0 | grep C19C | grep HPL")
    cat(fn,"\n")
    kin <- data.frame(fread(cmd=fn,header=F))
    kin <- kin %>% group_by(V2) %>% top_n(wt=V6,n=1)
    nn <- kin$V1
    names(nn) <- kin$V2
    nn[kin$V6<0] <- NA
    dd$assig2 <- nn[dd$assignment]
    dd  
},mc.cores=10)

demux <- do.call(rbind,demux)
###

### output
opfn <- paste0(outdir,"1_souporcell.ALL.rds")
write_rds(demux, opfn)


opfn <- paste0(outdir,"1_souporcell.SNG.rds")
demux %>% filter(status=="singlet") %>%
    write_rds(opfn)

table(demux$status)

table(demux$assig2)
