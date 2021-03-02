######################################
### cell counts ###
######################################


library(tidyverse)
library(parallel)
##library(data.table)

### 
outFolder="./1_souporcell_output/"

opfn <- "./1_souporcell_output/1_souporcell.ALL.rds"
demux <- read_rds(opfn)


aa <- demux %>% filter(!is.na(assig2),status=="singlet") %>%
    select(barcode,status,Sample_ID=assig2) %>%
    mutate(EXP=gsub("_.*","",barcode))


##SNG.BEST.GUESS
##
cc <- read_tsv("/nfs/rprdata/scilab/labor2/Covid19.Samples.txt") %>%
    mutate(Sample_ID=paste(Pregnancy_ID,Origin,sep="-"))
head(cc)


aa <- aa %>% select(EXP,Sample_ID) %>% left_join(cc) %>%
    select(Sample_ID,Pregnancy_ID,Origin,Condition,FetalSex,EXP)


cell.counts <- aa %>% group_by(EXP,Sample_ID,Pregnancy_ID,Origin,Condition) %>%
    summarize(n=n()) %>%
    filter(n>50)

write_tsv(cell.counts,paste0(outFolder,"cell.counts.tsv"))
