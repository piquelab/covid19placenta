
rm(list=ls())
#library(hta20probeset.db)
library(limma)
library(annaffy)
library(annotate)
library(lme4)
library(marray)
library(splines)
library(MASS)
library(ROCR)
require(pROC)
library(epiR)
library(org.Hs.eg.db)


cell.type.annotation<-read_tsv("cell.type.annotation.v3.tsv")
clust2Names<-cell.type.annotation$Potential.final #c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Stromal Fibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(cell.type.annotation$Cluster,":",clust2Names)
names(clust2Names)<-as.character(cell.type.annotation$Cluster)
clust2Names["4"]<-"4:Macrophage-2 (Hofbauer)" 

# load("GSE96083/eset_gse96083.RData")
# eset_gse96083
# load("GSE96083/ano_GSE96083.RData")

  
load("GSE96083/anoSC2_eset_GSE96083.RData")  
#outFolder<-"./14_GSE96083_analysis/"
outFolder<-"./14_GSE96083_DEGsig_analysis/"
system(paste0("mkdir -p ",outFolder))


#load("anoSC2_eset_GSE96083.RData")
ano=anoSC2
eset=as.matrix(esetSC2)

#load("eset_SC3_v11.RData")
#ano=anoSC3
#eset=as.matrix(esetSC3)
ano$Int=factor(cut(ano$GA,breaks=c(23,34)))


anpack="org.Hs.eg.db"
rownames(eset)<-gsub("_at","",rownames(eset))

SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
eset=eset[!is.na(SYMBOLS),]

SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
eset=eset[order(apply(eset,1,mean),decreasing=TRUE),]
SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
eset=eset[!duplicated(SYMBOLS),]
SYMBOLS<-unlist((lookUp(rownames(eset), anpack, 'SYMBOL')))
rownames(eset)<-SYMBOLS


#SCGeneSets=read.table("MarkersAll.txt",sep=" ",header=TRUE,stringsAsFactors = FALSE)



m2 = read_tsv("5_harmonyClusters_withcovid19control_DGE/ClusterDEG.tsv")
eg = bitr(m2$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="symbol"
head(eg)
# non-na ENTREZID were included
m2 <- m2 %>% left_join(eg) %>% filter(!is.na(ENTREZID))
m3 <- m2 %>% filter(p_val_adj<0.1,avg_log2FC>0.5) %>%
  group_by(gene) %>%
  mutate(H=log2(length(cluster))) %>%
  filter(H<=1) %>%
  ungroup()

table(m3$cluster)
dat <- m3 %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC) %>% ungroup()
table(dat$cluster)
dat$cluster<-clust2Names[as.character(dat$cluster)]




##### DEGs as signatures

res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_with_covidcontrol_res1.0_library/ALL.combined.2022-03-29.tsv")
# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","cluster","Origin"),sep="_",remove=FALSE) #Cell_type
res<-res %>% filter(!is.na(padj) & padj<0.1)
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
dat <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

# up-regulated DEGs
#dat <- dat %>% group_by(cluster) %>% top_n(n = 20, wt = log2FoldChange) %>% ungroup()

# both up-regulated and down-regulated DEGs
dat <- dat %>% group_by(cluster) %>% top_n(n = 20, wt = log2FoldChange) %>% ungroup()
table(dat$cluster)
dat$cluster<-clust2Names[as.character(dat$cluster)]
dat$cluster<-paste0(dat$cluster,"_",dat$Origin)
dat<-dat %>% dplyr::select(cluster,ENTREZID,symbol=gene_name,avg_log2FC=log2FoldChange,p_val=pvalue,p_val_adj=padj)


SCGeneSets<-dat

pile=NULL
for(a in unique(SCGeneSets$cluster)){
  tmp=SCGeneSets[SCGeneSets$cluster==a,]
  tmp=tmp[tmp$symbol%in%rownames(eset),]
  tmp=tmp[order(tmp$avg_log2FC,decreasing=TRUE),]
  pile=rbind(pile,tmp[1:min(c(25,dim(tmp)[1])),c("symbol","cluster")])
}
names(pile)<-c("Symbol","Type")
genes=pile


  
  eset2=NULL
  for(mi in unique(ano$Int)){
    esetM=apply(eset[,ano$Group=="Control"&ano$Int==mi],1,mean)
    esetS=apply(eset[,ano$Group=="Control"&ano$Int==mi],1,sd)
    eset2=cbind(eset2,apply(eset[,ano$Int==mi],2,function(x){(x-esetM)/esetS}))
    #eset2=cbind(eset2,apply(esett[,ano$Int==mi],2,function(x){(x-esetM)}))
  }
  
  #eset2=esett
  
  
  celsig=list()
  for(cs in unique(genes$Type)){
    celsig[[cs]]<-intersect(genes$Symbol[genes$Type==cs],SYMBOLS)
  }
  
  toad=NULL
  rw<-c()
  for (cs in names(celsig)){
    print (cs)
    
    if (any(SYMBOLS%in%celsig[[cs]]))
    {
      toad=rbind(toad,apply(as.matrix(eset2[SYMBOLS%in%celsig[[cs]],]),2,mean))
      rw<-c(rw,cs)
    }
      
      
    }
  rownames(toad)<-rw
  
  mat=toad
  colnames(mat)<-ano$Group
  # htm=heatmap(mat,Colv=NA,col=maPalette(low = "yellow", high = "blue", k = 50),
  #             ColSideColors=ifelse(ano$Group=="sPTD","orange","lightblue"))
  # library(Heatplus)
  # heatmap_2(mat, legend=1, legfrac=7,Colv=NA,col=maPalette(low = "yellow", high = "blue", k = 50))
  # 
  
  #save(esett,eset2, toad,ano,celsig,SYMBOLS,file="forpredPTL.RData")
  

  
  tg=data.frame(ID=rownames(toad),ID2=rownames(toad))
  
xx=toad
nm=paste0(outFolder,"PRB_sc-sigs_sPTD_23-34_revised.pdf")

    nint=length(levels(ano$Int))
    pdf(nm)
    par(bg="white")
    psall=NULL
    
    for(an in as.character(tg$ID)){
      ano2=ano
      ano2$Y=as.vector(as.matrix((xx[an,ano2$SampleID])))
      rg=range(ano2$Y);rg[1]<-rg[1]-0.2;rg[2]<-rg[2]+0.5
         
      levs=levels(ano2$Int)
      ps=NULL
      ndiag=NULL
      for(i in 1:nint){
        ps=c(ps,wilcox.test(Y~Group,data=ano2[ano2$Int%in%c(levs[i]),])$p.value)
        #ps=c(ps,t.test(Y~Group,data=ano2[ano2$Int%in%c(levs[i]),],var.equal =TRUE)$p.value)
              #}
      }
      
      
      psall=c(psall,ps)
      print(ps)
      #if(ps<0.0099){ #0.00765
      if(ps<0.2){ #0.00765
        boxplot(Y~Int,ano2[ano2$Group=="Control",],boxwex=0.25,col="lightblue",main=an,ylim=rg,
                xlab="Gestational age (weeks)",ylab="Average Z-score",cex.lab=1.3)
        
        boxplot(Y~Int,ano2[ano2$Group=="sPTD",],boxwex=0.25,col="orange",add=TRUE,at=(1:nint)+0.3,axes=FALSE)
        
      text(1:nint,y=rep(max(ano2$Y+0.1),nint),labels=paste("p=",round(ps,4),sep=""),cex=0.7)
      text(1:nint,y=rep(min(ano2$Y)-0.15,nint),labels=paste("n1=",table(ano$Int[ano$Group=="Control"]),sep=""),cex=0.6)
      text(1:nint+0.3,y=rep(min(ano2$Y)-0.15,nint),labels=paste("  n2=",table(ano$Int[ano$Group!="Control"]),sep=""),cex=0.6)
        #legend("top",horiz=TRUE,cex=0.8,fill=c("lightblue","orange",NA),pch=c(NA,NA,17),col=c(NA,NA,"red"),legend=c("Normal Pregnancy","Early PE","Post Diagnosis"))
      legend("top",horiz=TRUE,cex=0.8,fill=c("lightblue","orange"),legend=c("Normal Pregnancy","sPTD"))
      #abline(h=0)
      }
      #plot(Y~GA,ano2[ano2$Group=="Control",])
      #points(Y~GA,ano2[ano2$Group=="sPTD",],col="red")
    }
    dev.off()
    
 psa=p.adjust(psall,"fdr")
tg$ID[psa<0.1]
