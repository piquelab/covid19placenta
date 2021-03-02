###################################################
### pathway / gene ontology enrichment analysis ###
###################################################


library(tidyverse)
library(qqman)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(stringr)

# filtering 
filter_sample<-c("HPL20874") #inflamation
outFolder <- paste0("./11_pathway_enrichment_",paste(filter_sample,collapse="_"),"/")
system(paste0("mkdir -p ",outFolder))


#HPL20874
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_filtered_HPL20874_2020-11-28/ALL.combined.2020-11-28.tsv")

# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)

# Removing na pvalues
# Grouping pvalues based on the Location,Cell_type,and Origin
# Adding a column showing the rank of each pvalue devided by the number of pvalues in each group 
res2 <- res %>% filter(!is.na(pvalue)) %>%
    arrange(pvalue) %>%
    group_by(Location,Cell_type,Origin) %>%
    mutate(r=rank(pvalue, ties.method = "random"),pexp=r/length(pvalue))

#cluster colors
new_names <- read_tsv("../2020-10-02/5_harmony_cellClass_plots_res0.8/ClusterAssignment.res0.8.tsv")
clust2Names <- new_names$scLabor_ID
names(clust2Names) <- new_names$seurat_clusters
cc <- new_names %>% select(scLabor_ID,color) %>% unique 
cluster.Colors <- cc$color
names(cluster.Colors) <- cc$scLabor_ID


#ENTREZID id 
eg = bitr(res2$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID


# non-na ENTREZID were included
res3 <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))


#exploratory analysis
# # counting the number of DE genes per cname   
DE_per_cname<-sapply(unique(res3$cname), function(x,padj_cutoff=0.1,log2FoldChange_cutoff=0.25){
    aux <- res3 %>% filter(cname==x)
    genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    length(genes)
    
})

which(DE_per_cname>0)
cname_selected<-names(DE_per_cname)[which(DE_per_cname>0)]


#qvalueCutoff  = 0.05

pathway_enrich<-function(res_gene=res3,cname_select="CAM_T-cell_M",padj_cutoff=0.1,log2FoldChange_cutoff=0.25)
{
    pathway_enrich_cname_dir<-paste0(outFolder,cname_select,"/")
    system(paste0("mkdir -p ",pathway_enrich_cname_dir))
    result<-list()
    aux <- res_gene %>% filter(cname==cname_select)
    genes <- filter(aux,padj<padj_cutoff,abs(log2FoldChange)>log2FoldChange_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
    geneUniv <- aux %>% dplyr::select(ENTREZID) %>% unlist
    ##geneList <- aux$log2FoldChange
    geneList <- -log10(aux$pvalue)
    names(geneList) <- aux$ENTREZID
    geneList = sort(geneList, decreasing = TRUE)
    message(".................................")
    message("Number of DE genes: ",length(genes))
    #print(length(genes))
    
    message(".................................")
    message("enrichGO")
    ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
    print(head(ego))
    result$enrichGO<-ego
    save(ego,file=paste0(pathway_enrich_cname_dir,"ego.RData"))
    write.csv(ego,file=paste0(pathway_enrich_cname_dir,"ego.csv"))
    
    print(".................................")
    print("enrichKEGG")
    ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
    print(head(ekegg)) 
    result$enrichKEGG<-ekegg
    save(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.RData"))
    write.csv(ekegg,file=paste0(pathway_enrich_cname_dir,"ekegg.csv"))
    
    
    message(".................................")
    message("enrichPathway")
    erpath <- enrichPathway(gene=genes,universe=geneUniv)
    print(head(erpath))
    result$enrichPathway<-erpath
    save(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.RData"))
    write.csv(erpath,file=paste0(pathway_enrich_cname_dir,"erpath.csv"))
    
    message(".................................")
    message("gseGO")
    # BP: biological_process, CC: cellular_component, MF: molecular_function
    gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP")
    print(head(gseGO.res))
    result$gseGO<-gseGO.res
    save(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.RData"))
    write.csv(gseGO.res,file=paste0(pathway_enrich_cname_dir,"gseGO.res.csv"))
    
    
    message(".................................")
    message("gsePathway")
    gseRPath.res <- gsePathway(geneList)
    print(head(gseRPath.res))
    result$gsePathway<-gseRPath.res
    save(gseRPath.res,file=paste0(pathway_enrich_cname_dir,"gseRPath.res.RData"))
    write.csv(gseRPath.res,file=paste0(pathway_enrich_cname_dir,"gseRPath.res.csv"))
    return (result)
    }


result_pathway_en_list<-lapply(cname_selected, function(x) return(pathway_enrich(res3,x)))
names(result_pathway_en_list)<-cname_selected

save(result_pathway_en_list,file=paste0(outFolder,"pathwayEnrich_result.RData"))
which(DE_per_cname>0)


##########################################################################################  
#####                                  dot plot
##########################################################################################

#per cname

# how to calculate GeneRatio
#https://github.com/YuLab-SMU/DOSE/issues/20

filter_sample<-c("HPL20874") #inflamation
outFolder <- paste0("./11_pathway_enrichment_",paste(filter_sample,collapse="_"),"/")

system(paste0("mkdir -p ",outFolder))

load(paste0(outFolder,"pathwayEnrich_result.RData"))
cname_selected<-names(result_pathway_en_list)

outFolder_cname_plots <- paste0("./12_pathway_enrichment_",paste(filter_sample,collapse="_"),"_plots/")
system(paste0("mkdir -p ",outFolder_cname_plots))

##selecting cname with #DE >5
DE_per_cname_select<-names(DE_per_cname)[which(DE_per_cname>=5)]
result_pathway_en_list<-result_pathway_en_list [DE_per_cname_select]
cname_selected<-names(result_pathway_en_list)
#gseGO
res_gseGO_list<-lapply(cname_selected, function(x)
{
    
    rs<-result_pathway_en_list[[x]]$gseGO@result %>% filter(qvalues<=0.05)
    dim1<-dim(rs)[1]

    if(min(dim1,5)>0)
    {
        res_en<-rs
        
        # to calculate GeneRatio=count/setSize
        
        #count
        gene_count<- res_en %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
        
        ## merge with the original dataframe
        dot_df<- left_join(res_en, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
        
        dot_df<-dot_df[1:min(dim1,5),c("ID","Description" ,"enrichmentScore","p.adjust","GeneRatio")]
        dot_df$cname<-rep(x,min(dim1,5))
        dot_df
        }
    })
    

res_df_gseGO <- do.call(rbind,res_gseGO_list)

# Count is the number of core genes and GeneRatio is Count/setSize

pdf(paste0(outFolder_cname_plots,"gseGO_cname_DotPlot.pdf"),width=15,height=15)
ggplot(res_df_gseGO, 
       aes(x = cname, y = Description)) + 
    geom_point(aes(size = enrichmentScore, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.03))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term" #enrichmentScore
    ylab(NULL)+ 
    xlab(NULL)
dev.off()


#enrichGO
res_enrichGO_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichGO)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichGO@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,10)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,10))
        }
        res_en }
}
    )  
    
res_df_enrichGO <- do.call(rbind,res_enrichGO_list)




####################

#enrichKEGG
res_enrichKEGG_list<-lapply(cname_selected, function(x)
{
    if(length(result_pathway_en_list[[x]]$enrichKEGG)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichKEGG@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,10)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,10))
        }
        res_en }
}
)  

res_df_enrichKEGG <- do.call(rbind,res_enrichKEGG_list)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

pdf(paste0(outFolder_cname_plots,"enrichKEGG_cname_DotPlot.pdf"),width=10,height=15)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = cname, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.1))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)
dev.off()


####################

#enrichPathway

res_enrichPathway_list<-lapply(cname_selected, function(x)
{
    
    
    if(length(result_pathway_en_list[[x]]$enrichPathway)>0)
    {
        rs<-result_pathway_en_list[[x]]$enrichPathway@result %>% filter(qvalue<=0.05)
        
        dim1<-dim(rs)[1]
        res_en<-NULL
        if(min(dim1,10)>0)
        {
            res_en<-rs
            res_en<-res_en[1:min(dim1,10),c("ID","Description" ,"GeneRatio","p.adjust")]
            res_en$cname<-rep(x,min(dim1,10))
        }
        res_en }
}
)  

res_df_enrichPathway <- do.call(rbind,res_enrichPathway_list)

res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})

############################################################
### Combined DEGs 
############################################################

genes <- filter(res3,padj<0.1,abs(log2FoldChange)>0.25) %>% dplyr::select(ENTREZID) %>% unlist %>% unique
geneUniv <- res3 %>% dplyr::select(ENTREZID) %>% unlist %>% unique

geneList <- -log10(res3$pvalue)
names(geneList) <- res3$ENTREZID
geneList = sort(geneList, decreasing = TRUE)


message(".................................")
message("enrichGO")
ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
print(head(ego))


print(".................................")
print("enrichKEGG")
ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
print(head(ekegg)) 

message(".................................")
message("enrichPathway")
erpath <- enrichPathway(gene=genes,universe=geneUniv)
print(head(erpath))


message("gsePathway")
gseRPath.res <- gsePathway(geneList)
print(head(gseRPath.res))

save(gseRPath.res,file=paste0(outFolder,"gsePathway.res_combined.RData"))


res_df<-gseRPath.res@result[1:10,]
res_df<-res_df %>% filter(qvalues<=0.05)
pdf(paste0(outFolder_cname_plots,"gsePathway.res_combined_DotPlot.pdf"),width=20,height=10)
ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
       aes(x = enrichmentScore, y = Description)) + 
    geom_point(aes(size = enrichmentScore, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.03))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
    ylab(NULL) 
dev.off()


message(".................................")
message("gseGO")
gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP")
print(head(gseGO.res))
save(gseGO.res,file=paste0(outFolder,"gseGO.res_combined.RData"))

#plot
res_df<-gseGO.res@result[1:10,]
res_df<-res_df %>% filter(qvalues<=0.05)
pdf(paste0(outFolder_cname_plots,"gseGO.res_combined_DotPlot.pdf"),width=20,height=10)
ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
       aes(x = enrichmentScore, y = Description)) + 
    geom_point(aes(size = enrichmentScore, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.03))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
    ylab(NULL) 
dev.off()


message(".................................")
message("gsePathway")
gseRPath.res <- gsePathway(geneList)
print(head(gseRPath.res))
save(gseRPath.res,file=paste0(outFolder,"gseRPath.res_combined.RData"))


#plot
res_df<-gseRPath.res@result[1:10,]
res_df<-res_df %>% filter(qvalues<=0.05)
pdf(paste0(outFolder_cname_plots,"gseRPath.res_combined_DotPlot.pdf"),width=30,height=20)
ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
       aes(x = enrichmentScore, y = Description)) + 
    geom_point(aes(size = enrichmentScore, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.03))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
    ylab(NULL) 
dev.off()

########################################################
#comparison between PVBP_T-cell_M and CAM_T-cell_M
########################################################

filter_sample<-"HPL20874"
#loading all DE genes
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_filtered_HPL20874_2020-11-28/ALL.combined.2020-11-28.tsv")

eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)
e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID
res <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID))

# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))

res_PVBP_Tcell_M <- res %>% filter(cname=="PVBP_T-cell_M")

res_CAM_Tcell_M <- res %>% filter(cname=="CAM_T-cell_M")


# DE genes 

pvalue_cutoff=0.01
log2FoldChange_cutoff=0.1

CAM_Tcell_M_de<-filter(res_CAM_Tcell_M,pvalue<pvalue_cutoff) %>% dplyr::select(ENTREZID) %>% unlist
PVBP_Tcell_M_de<-filter(res_PVBP_Tcell_M,pvalue<pvalue_cutoff) %>% dplyr::select(ENTREZID) %>% unlist

length(PVBP_Tcell_M_de)
length(CAM_Tcell_M_de)

genes <- intersect(CAM_Tcell_M_de,PVBP_Tcell_M_de)  
geneUniv <- res$ENTREZID  

message(".................................")
message("Number of DE genes: ",length(genes))
#print(length(genes))

message(".................................")
message("enrichGO")
ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
print(head(ego))


print(".................................")
print("enrichKEGG")
ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
print(head(ekegg)) 
result$enrichKEGG<-ekegg


message(".................................")
message("enrichPathway")
erpath <- enrichPathway(gene=genes,universe=geneUniv)
print(head(erpath))

####################################################################################
# based on the shared DE genes between T-cell and reference
####################################################################################

filter_sample<-"HPL20874"
#loading all DE genes
res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_filtered_HPL20874_2020-11-28/ALL.combined.2020-11-28.tsv")

outFolder_cname_plots <- paste0("./12_pathway_enrichment_",paste(filter_sample,collapse="_"),"_plots/")
system(paste0("mkdir -p ",outFolder_cname_plots))


# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))

#ENTREZID id 
eg = bitr(res$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(eg)[1]="gene_name"
head(eg)

e2g <- eg$gene_name
names(e2g) <- eg$ENTREZID

# non-na ENTREZID were included
res3 <- res %>% left_join(eg) %>% filter(!is.na(ENTREZID) & !is.na(pvalue))
res_CAM_Tcell_M <- res3 %>% filter(cname=="CAM_T-cell_M")

Tref <- read_tsv("/nfs/rprdata/scilab/labor2/covid-analysis/TcellRef/TcellCovidCluster5.txt")
colnames(Tref)[1:4] <- c("gene_name","R.Log2FC","Rpval","Rpadj")
Tref <- Tref %>% filter(!is.na(R.Log2FC))

Tref <- Tref %>% left_join(eg) %>% filter(!is.na(ENTREZID))



res_CAM_Tcell_M <- res3 %>% filter(cname=="CAM_T-cell_M")

##########################################
#overlap between CAM and ref (based on padj<0.1,abs(log2FoldChange)>0.0 )
##########################################
res_CAM_Tcell_M_DE_adjust<-res_CAM_Tcell_M %>% filter(padj<0.1,abs(log2FoldChange)>0.0) 
Tref_DE_adjust<-Tref %>% filter(Rpadj<0.1,abs(R.Log2FC)>0.0) 
de_adjust_cam_not_ref<-res_CAM_Tcell_M_DE_adjust$ENTREZID [which(!res_CAM_Tcell_M_DE_adjust$ENTREZID %in% Tref_DE_adjust$ENTREZID)]
length(de_adjust_cam_not_ref)

de_adjust_cam_not_ref_sig<-res_CAM_Tcell_M_DE_adjust %>% filter(!ENTREZID %in% Tref_DE_adjust$ENTREZID)

#
write.csv(de_adjust_cam_not_ref_sig,file = paste0(outFolder,"de_cam_not_ref.csv"))

##########################################
#overlap between DE genes from CAM and refrence based on pvalue 
##########################################
res_CAM_Tcell_M_DE<-res_CAM_Tcell_M %>% filter(pvalue<0.01) 
Tref_DE<-Tref %>% filter(Rpval<0.01) 
Tref_DE_gene<-Tref_DE$ENTREZID
res_CAM_Tcell_M_DE_gene<-res_CAM_Tcell_M_DE$ENTREZID
overlap_DE_gene<-intersect(res_CAM_Tcell_M_DE_gene,Tref_DE_gene)
length(overlap_DE_gene)

#overlap DE genes and re-adjusting by "fdr"
res_DE_pvalue_overlap<-res_CAM_Tcell_M %>% filter(ENTREZID %in% overlap_DE_gene) 
res_DE_pvalue_overlap$padj<-p.adjust(res_DE_pvalue_overlap$pvalue ,"fdr")

genes <- res_DE_pvalue_overlap %>%filter(padj<0.1,abs(log2FoldChange)>0.25) %>% dplyr::select(ENTREZID) %>% unlist %>% unique
geneUniv <- res_CAM_Tcell_M %>% dplyr::select(ENTREZID) %>% unlist %>% unique

geneList<-res_CAM_Tcell_M$pvalue
geneList<--log10(geneList)
names(geneList)<-res_CAM_Tcell_M$ENTREZID #%%in% overlap_DE_gene
geneList = sort(geneList, decreasing = TRUE)

message(".................................")
message("enrichGO")
ego <- enrichGO(gene=genes,universe=geneUniv, OrgDb=org.Hs.eg.db,ont="BP")
print(head(ego))
###


res_df_enrichGO<-ego@result %>% filter(qvalue<=0.05)

res_df_enrichGO$GeneRatio<-sapply(res_df_enrichGO$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})
pdf(paste0(outFolder_cname_plots,"enrichGO.res_CAM_Tcell_M_overlapDE_padjust_DotPlot.pdf"),width=15,height=10)
ggplot(res_df_enrichGO, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.032))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)+
    xlab("GeneRatio") 
    #coord_fixed(ratio = 1)
#ggtitle("GO pathway enrichment")
dev.off()

print(".................................")
print("enrichKEGG")
ekegg <- enrichKEGG(gene=genes,universe=geneUniv,organism="hsa")
print(head(ekegg)) 

res_df_enrichKEGG<-ekegg@result%>% filter(qvalue<=0.05)

res_df_enrichKEGG$GeneRatio<-sapply(res_df_enrichKEGG$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})
pdf(paste0(outFolder_cname_plots,"enrichKEGG.res_CAM_Tcell_M_overlapDE_padjust_DotPlot.pdf"),width=5,height=5)
ggplot(res_df_enrichKEGG, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.03))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=10)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)
dev.off()


message(".................................")
message("enrichPathway")
erpath <- enrichPathway(gene=genes,universe=geneUniv)
print(head(erpath))

res_df_enrichPathway<-erpath@result%>% filter(qvalue<=0.05)

res_df_enrichPathway$GeneRatio<-sapply(res_df_enrichPathway$GeneRatio, function(x){
    numden<-unlist(strsplit(x,"/"))
    return (as.numeric(numden[1])/as.numeric(numden[2]))
})
pdf(paste0(outFolder_cname_plots,"enrichPathway.res_CAM_Tcell_M_overlapDE_padjust_DotPlot.pdf"),width=10,height=5)
ggplot(res_df_enrichPathway, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Description)) + 
    geom_point(aes(size = GeneRatio, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.00001, 0.03))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=20)) +
    labs(size="GeneRatio",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)+ 
    xlab(NULL)
dev.off()



message(".................................")
message("gseGO")
gseGO.res <- gseGO(geneList,  OrgDb=org.Hs.eg.db,ont="BP")
print(head(gseGO.res))

res_df_gseGO<-gseGO.res@result %>% filter(qvalues<=0.05)

#plot
res_df<-res_df_gseGO

pdf(paste0(outFolder_cname_plots,"gseGO.res_CAM_Tcell_M_overlapDE_padjust_DotPlot.pdf"),width=22,height=12)
ggplot(res_df, # you can replace the numbers to the row number of pathway of your interest
       aes(x = enrichmentScore, y = Description)) +
    geom_point(aes(size = enrichmentScore, color = p.adjust)) +
    theme_bw(base_size = 11) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradient(low = "red",  high = "blue", space = "Lab", limit = c(0.001, 0.05))+
    theme(axis.text.x = element_text(angle = 45,hjust=1),text = element_text(size=30)) +
    labs(size="enrichmentScore",color="p.adjust") + #x="",y="GO term"
    ylab(NULL)
dev.off()

message(".................................")
message("gsePathway")
gseRPath.res <- gsePathway(geneList)
print(head(gseRPath.res))