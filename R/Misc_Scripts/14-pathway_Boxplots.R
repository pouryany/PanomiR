rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyverse)
options(stringsAsFactors = F)

covariates <- read.csv('../Data/TCGA-LIHC-COV.csv', row.names = 1)
covariates <- covariates[,c("barcode","shortLetterCode")]
condition = 'shortLetterCode'
out.dir   = '../test_cases/LIHC/Output/'
data.dir  = '../test_cases/LIHC/Data/'



de.paths0 <- readRDS(paste0(data.dir,"Pathways_object.RDS"))
de.paths  <- de.paths0$DEP

top.paths <- rownames(de.paths[1:9,])

path.residuals <- de.paths0$PathwayResiduals
path.residuals <- path.residuals[top.paths,]
path.residuals <- t(path.residuals)

library(reshape2)

path.residuals          <- reshape2::melt(path.residuals)
names(path.residuals)   <- c("barcode","Pathway","Expression")
path.residuals          <- left_join(path.residuals,covariates)



path.residuals$Pathway  <- gsub("Pathway.","",path.residuals$Pathway)
path.residuals$Pathway  <- gsub("^([A-Z]*)_","\\1: ",path.residuals$Pathway)
path.residuals$Pathway  <- gsub("_"," ",path.residuals$Pathway)



library(cowplot)
ggplot(data = path.residuals, aes(x=shortLetterCode, y=Expression)) + 
    geom_boxplot(aes(fill=shortLetterCode)) +
    facet_wrap( ~ Pathway, scales="free") + 
    labs(fill = "Tissue type") + 
    theme_bw()

de.paths <- de.paths0$DEP


de.paths <- de.paths %>% rownames_to_column(.,var = "Pathway") %>% mutate_if(.,is.numeric, signif, digits = 3)


de.paths = de.paths %>% 
    dplyr::mutate(Direction = logFC/abs(logFC),
                  Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                  Direction = as.character(Direction))
de.paths[de.paths$adj.P.Val > 0.05,]$Direction = 'NONE'



de.paths$Pathway_Name  <- gsub("Pathway.","",de.paths$Pathway)
de.paths$Pathway_Name  <- gsub("^([A-Z]*)_","\\1: ",de.paths$Pathway_Name)
de.paths$Pathway_Name  <- gsub("_"," ",de.paths$Pathway_Name)

clustering_louvain        <- read.csv("../temp22/Pathways_cluster_louvain.csv",row.names = "X")
clustering_louvain        <- left_join(clustering_louvain,de.paths)
write.csv(de.paths,"../temp_figures/DE_pathways_rep_LIHC.csv")
write.csv(clustering_louvain,"../temp_figures/cluster_DE_pathways_rep_LIHC.csv")




enrichment_res <- read.csv("../test_cases/LIHC/Output/Enrichment_pathways.csv")
names(enrichment_res)  <- paste0("Enrichment_",names(enrichment_res))


names(enrichment_res)[1] <- "Pathway_Name"


enrichment_res <- left_join(de.paths,enrichment_res)

enrichment_res[is.na(enrichment_res$Enrichment_BgRatio),
               ]$Enrichment_BgRatio <- "0"

enrichment_res[is.na(enrichment_res$Enrichment_p.adjust),
               ]$Enrichment_p.adjust <- 1

enrichment_res[is.na(enrichment_res$Enrichment_Count),
               ]$Enrichment_Count <- 0


write.csv(enrichment_res,"../temp_figures/full_pathways_LIHC.csv")




all.res  <- list.files("../test_cases/LIHC/Output/cluster_louvain_Prioritization_Tarbase/",full.names = T)


for(res in all.res){
    res.tab0 <- read.csv(res,row.names = "X")
    
    res.tab0[is.na(res.tab0$cluster_hits),]$cluster_hits <- 0
    
    res.tab0 <- res.tab0 %>% mutate_if(., is.numeric,signif,digits = 3)
    
    res.new  <- gsub("../test_cases/LIHC/Output/cluster_louvain_Prioritization_Tarbase/","../temp_figures/",res)
    res.new  <- gsub("x2","louvain_tarbase",res.new)
    print(nrow(res.tab0[res.tab0$AggInv_fdr< 0.00001,]))
    write.csv(res.tab0,res.new)
    
}


all.res  <- list.files("../test_cases/LIHC/Output/cluster_louvain_Prioritization_TargetScan03/",full.names = T)


for(res in all.res){
    res.tab0 <- read.csv(res,row.names = "X")
    if(any(is.na(res.tab0$cluster_hits)))
    res.tab0[is.na(res.tab0$cluster_hits),]$cluster_hits <- 0
    
    res.tab0 <- res.tab0 %>% mutate_if(., is.numeric,signif,digits = 3)
    
    res.new  <- gsub("../test_cases/LIHC/Output/cluster_louvain_Prioritization_TargetScan03/","../temp_figures/",res)
    res.new  <- gsub("x2","louvain_targetScan_03",res.new)
    
    write.csv(res.tab0,res.new)
    
}









