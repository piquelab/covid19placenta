###################################################
### QQ plot                                     ###
###################################################

library(tidyverse)
##library(knitr)
library(DESeq2)
##library(annotables)
library(qqman)

filter_sample<-c("HPL20874")
outFolder <- paste0("./8_outputs_DESeq_Plots_",paste(filter_sample,collapse="_"),"/")
#outFolder <- paste0("./8_outputs_DESeq_Plots/")
system(paste0("mkdir -p ",outFolder))


res <- read_tsv("./7_outputs_DESeq_ConditionsByCluster_filtered_HPL20874_HPL20875_2020-11-27/ALL.combined.2020-11-28.tsv")

## To do next: Plot top genes or heatmap, sumarize effect sizes w/ forest plot.
## pathway analysis, enrichment analysis. 

fname=paste0(outFolder,"all.qqplot.png");
png(fname,width=800,height=800)
qq(res$pvalue)
dev.off()

# Adding location, cell type, and origin columns 
res <- res %>% separate(cname,c("Location","Cell_type","Origin"),sep="_",remove=FALSE)


###################################################
# Grouping pvalues based on the Location,Cell_type,and Origin
# Adding a column showing the rank of each pvalue devided by the number of pvalues in each group 
###################################################
res2 <- res %>% filter(!is.na(pvalue)) %>%
    arrange(pvalue) %>%
    group_by(Location,Cell_type,Origin) %>%
    mutate(r=rank(pvalue, ties.method = "random"),pexp=r/length(pvalue))


new_names <- read_tsv("../2020-10-02/5_harmony_cellClass_plots_res0.8/ClusterAssignment.res0.8.tsv")
clust2Names <- new_names$scLabor_ID
names(clust2Names) <- new_names$seurat_clusters
cc <- new_names %>% select(scLabor_ID,color) %>% unique 
cluster.Colors <- cc$color
names(cluster.Colors) <- cc$scLabor_ID
tempcol<-cluster.Colors["B-cell"]
cluster.Colors["B-cell"]<-cluster.Colors["T-cell"]
cluster.Colors["T-cell"]<-tempcol

###################################################
# qqplot to show the p-values splited by Origin and Location  
###################################################
fname=paste0(outFolder,"split.qqplot.png");
p1 <- res2 %>%
    ggplot(aes(x=-log10(pexp),y=-log10(pvalue),color=Cell_type)) +
    geom_point() +
    scale_color_manual(values=cluster.Colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell Type")) +
    geom_abline(slope=1,intercept=0) +
    facet_grid(Origin ~ Location) +
    xlab(expression(Expected -log[10](p))) +
    ylab(expression(Observed -log[10](p))) + 
    theme_bw()

ggsave(fname,p1,width=6,height=4.5)


