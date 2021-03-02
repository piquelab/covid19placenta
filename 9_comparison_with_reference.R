###################################################
### Comparison T-cell with reference            ###
###################################################
# reference: Meckiff, Benjamin J., et al. "Imbalance of regulatory and cytotoxic SARS-CoV-2-reactive CD4+ T cells in COVID-19." Cell 183.5 (2020): 1340-1353.

library(tidyverse)
library(DESeq2)
##library(annotables)
library(qqman)

### Compared to cell T-cell paper

filter_sample<-"HPL20874"
#loading all DE genes
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_filtered_HPL20874_2020-11-28/ALL.combined.2020-11-28.tsv")

outFolder<-"./8_outputs_DESeq_Plots/"
if(!is.null(filter_sample))
    outFolder <- paste0("./8_outputs_DESeq_Plots_",paste(filter_sample,collapse="_"),"/")

# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
#res <- res %>% filter(!is.na(pvalue))

###################################################
#loading reference 
###################################################
Tref <- read_tsv("/nfs/rprdata/scilab/labor2/covid-analysis/TcellRef/TcellCovidCluster5.txt")
colnames(Tref)[1:4] <- c("gene_name","R.Log2FC","Rpval","Rpadj")

Tref <- Tref %>% filter(!is.na(R.Log2FC))


res4 <- res %>% filter(cname=="CAM_T-cell_M")
dim(res4)

res_CAM_Tcell_M<-res4
resJoin <- res4 %>% inner_join(Tref) 

###################################################
# comparison with CAM_T-cell_M using all genes 
###################################################

cor.test(resJoin$log2FoldChange,resJoin$R.Log2FC,method="spearman",na.rm=TRUE)
cor.test(resJoin$log2FoldChange,resJoin$R.Log2FC,method="pearson",na.rm=TRUE)


###################################################
# comparison with CAM_T-cell_M based on common DE genes between CAM and reference 
###################################################
res5 <- res4 %>% filter(pvalue<0.01)
Tref2 <- Tref %>% filter(Rpval<0.01)
resJoin <- res5 %>% inner_join(Tref2) 
cor.test(resJoin$log2FoldChange,resJoin$R.Log2FC,method="spearman",na.rm=TRUE)

cor.test(resJoin$log2FoldChange,resJoin$R.Log2FC,method="pearson",na.rm=TRUE)
dim(resJoin)


p2 <- resJoin %>% arrange(-pvalue) %>% 
    ggplot(aes(R.Log2FC,log2FoldChange,color=pvalue<0.01)) +
    geom_point() +
    scale_color_manual(values=c("gray","black")) +
    theme_bw()

fname=paste0(outFolder,"Tcell.CAM.Val.png");
ggsave(fname,p2,width=6,height=4.5)


###################################################
#comparison with PVBP_T-cell_M
###################################################
res4 <- res %>% filter(cname=="PVBP_T-cell_M")
dim(res4)

res_PVBP_Tcell_M<-res4
resJoin <- res4 %>% inner_join(Tref) 

cor.test(resJoin$log2FoldChange,resJoin$R.Log2FC,method="spearman",na.rm=TRUE)

cor.test(resJoin$log2FoldChange,resJoin$R.Log2FC,method="pearson",na.rm=TRUE)


###################################################
# comparison with PVBP_T-cell_M and reference T-cell using common DE genes 
###################################################
res6 <- res4 %>% filter(pvalue<0.01)
Tref3 <- Tref %>% filter(Rpval<0.01)
resJoin <- res6 %>% inner_join(Tref3) 
cor.test(resJoin$log2FoldChange,resJoin$R.Log2FC,method="spearman",na.rm=TRUE)
cor.test(resJoin$log2FoldChange,resJoin$R.Log2FC,method="pearson",na.rm=TRUE)
dim(resJoin)


p2 <- resJoin %>% arrange(-pvalue) %>% 
    ggplot(aes(R.Log2FC,log2FoldChange,color=pvalue<0.01)) +
    geom_point() +
    scale_color_manual(values=c("gray","black")) +
    theme_bw()

fname=paste0(outFolder,"Tcell.PVBP.Val.png");
ggsave(fname,p2,width=6,height=4.5)


###################################################
#comparison between PVBP_T-cell_M and CAM_T-cell_M
###################################################
# all genes

res_PVBP_Tcell_M_fc_pvalue <- res_PVBP_Tcell_M %>% select(gene_name,log2FoldChange,pvalue,padj)
res_CAM_Tcell_M_fc_pvalue <- res_CAM_Tcell_M %>% select(gene_name,log2FoldChange,pvalue,padj)
colnames(res_CAM_Tcell_M_fc_pvalue)<-c("gene_name","CAMlog2FoldChange","CAMpvalue","CAMpadj")

resJoin <- res_CAM_Tcell_M_fc_pvalue %>% inner_join(res_PVBP_Tcell_M_fc_pvalue) 
cor.test(resJoin$log2FoldChange,resJoin$CAMlog2FoldChange,method="spearman",na.rm=TRUE)


resJoin$col <- "grey"
resJoin[!is.na(resJoin$CAMpvalue) & resJoin$CAMpvalue <0.01 , "col"] <- "black"
resJoin[!is.na(resJoin$pvalue) & resJoin$pvalue < 0.01 , "col"] <- "black"

p2 <- resJoin %>% arrange(-CAMpvalue) %>% 
    ggplot(aes(log2FoldChange,CAMlog2FoldChange,color=col)) +
    geom_point() +
    scale_color_manual(name="Pvalue < 0.01 (CAM or PVBP)",values=c("black","gray"),labels=c("True", "False")) +
    #scale_color_identity()+
    xlab("Log2FoldChange (PVBP)") + ylab("Log2FoldChange (CAM)")
theme_bw()

fname=paste0(outFolder,"Tcell.PVBP.CAM.Val.png");
ggsave(fname,p2,width=6,height=4.5)
