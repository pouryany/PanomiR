# Compare prioritization results to DE miRNAs. Trying to attribute miRNAs to cluster.
rm(list = ls())
library(dplyr)
library(tibble)


EXP.mirs  <- list.files("../test_cases/LIHC/Output/cluster_louvain_Prioritization_Tarbase/",full.names = T)
Pred.mirs <- list.files("../test_cases/LIHC/Output/cluster_louvain_Prioritization_TargetScan03/",full.names = T)




de.targets <- list()
for(res in 1:length(EXP.mirs)){
    
    predicted      <- read.csv(Pred.mirs[res],row.names = "X")
    sig.pred.mirs  <- predicted[predicted$AggInv_fdr<0.01,]

    
    res.tab0 <- read.csv(EXP.mirs[res],row.names = "X")
    res.tab <- res.tab0[res.tab0$AggInv_fdr<0.01,]
    
    de.vec  <- res.tab$x
    de.vec  <- paste(de.vec, collapse = "|")
    
    
    
    
    sel.mirs <- grep(de.vec,sig.pred.mirs$x,value = T)
    
    name.tag <- paste0("Cluster_",res)
    
    de.targets[[name.tag]]<- sort(sel.mirs)
    
    
}


de.target.tab <- list()
for(i in 1:length(de.targets)){
    temp.tab <- data.frame("Cluster" = names(de.targets)[i], "miRNA" = de.targets[[i]])
    de.target.tab <- rbind(de.target.tab,temp.tab)

}


write.csv(de.target.tab,"../temp_figures/Common_TarBase_TargetScan")

