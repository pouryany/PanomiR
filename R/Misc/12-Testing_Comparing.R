# Compare prioritization results to DE miRNAs. Trying to attribute miRNAs to cluster.
rm(list = ls())
library(dplyr)
library(tibble)


all.res  <- list.files("../test_cases/LIHC/Output/cluster_louvain_Prioritization_Tarbase/",full.names = T)
all.mirs <- readRDS("../Data/TCGA-LIHC-DEmiRNAs.RDS")
de.mirs  <- all.mirs[all.mirs$adj.P.Val<0.01,]

de.mirs <- rownames_to_column(de.mirs)
de.vec  <- de.mirs$rowname
de.vec  <- paste(de.vec, collapse = "|")

all.mirs <- rownames_to_column(all.mirs)
all.vec  <- all.mirs$rowname
all.vec  <- paste(all.vec, collapse = "|")



de.targets <- list()
n.targets  <- list()

for(res in all.res){
    res.tab0 <- read.csv(res,row.names = "X")
    res.tab <- res.tab0[res.tab0$AggInv_fdr<0.001,]
    
    all.ids  <- grep(all.vec,res.tab0$x,value = T)
    des.ids  <- grep(de.vec,res.tab0$x,value = T)
    
    sel.mirs <- grep(de.vec,res.tab$x,value = T)
    
    name.tag  <- unlist(strsplit(res,"/"))
    name.tag  <- tail(name.tag,1)
    name.tag  <- gsub(".csv","",name.tag)
    clust_tag <- tail(unlist(strsplit(name.tag,"_")),1)
        
    
    if(length(sel.mirs)>0)
    de.targets[[name.tag]] <- data.frame(cluster = clust_tag ,
                                         miRNA = sort(sel.mirs))
    n.targets[[name.tag]]  <- nrow(res.tab)
    
}

de.targets1 <- Reduce(rbind, de.targets)

write.csv(de.targets1, "../test_cases/LIHC/Output/DE_prioritzed_Tarbase.csv")
write.csv(all.mirs, "../test_cases/LIHC/Output/DE_miRs.csv")


a1 <- setdiff(de.targets$x2_LIHCGene_Tarbase1000_samples_clustNo_1$miRNA,
        de.targets$x2_LIHCGene_Tarbase1000_samples_clustNo_2$miRNA)


a2 <- setdiff(de.targets$x2_LIHCGene_Tarbase1000_samples_clustNo_2$miRNA,
        de.targets$x2_LIHCGene_Tarbase1000_samples_clustNo_1$miRNA)

write.csv(a1, "../test_cases/LIHC/Output/DE_prioritzed_Tarbase_Unique_1.csv")
write.csv(a2, "../test_cases/LIHC/Output/DE_prioritzed_Tarbase_Unique_2.csv")



all.res  <- list.files("../test_cases/LIHC/Output/cluster_louvain_Prioritization_TargetScan03/",full.names = T)
all.mirs <- readRDS("../Data/TCGA-LIHC-DEmiRNAs.RDS")
de.mirs  <- all.mirs[all.mirs$adj.P.Val<0.01,]

de.mirs <- rownames_to_column(de.mirs)
de.vec  <- de.mirs$rowname
de.vec  <- paste(de.vec, collapse = "|")

all.mirs <- rownames_to_column(all.mirs)
all.vec  <- all.mirs$rowname
all.vec  <- paste(all.vec, collapse = "|")



de.targets <- list()
n.targets  <- list()

for(res in all.res){
    res.tab0 <- read.csv(res,row.names = "X")
    res.tab <- res.tab0[res.tab0$AggInv_fdr<0.001,]
    
    all.ids  <- grep(all.vec,res.tab0$x,value = T)
    des.ids  <- grep(de.vec,res.tab0$x,value = T)
    
    sel.mirs <- grep(de.vec,res.tab$x,value = T)
    
    name.tag  <- unlist(strsplit(res,"/"))
    name.tag  <- tail(name.tag,1)
    name.tag  <- gsub(".csv","",name.tag)
    clust_tag <- tail(unlist(strsplit(name.tag,"_")),1)
    
    
    if(length(sel.mirs)>0)
        de.targets[[name.tag]] <- data.frame(cluster = clust_tag ,
                                             miRNA = sort(sel.mirs))
    n.targets[[name.tag]]  <- nrow(res.tab)
    
}

a1 <- setdiff(de.targets$x2_LIHCGene_TargetScan031000_samples_clustNo_1$miRNA,
              union(de.targets$x2_LIHCGene_TargetScan031000_samples_clustNo_2$miRNA,
                    de.targets$x2_LIHCGene_TargetScan031000_samples_clustNo_3$miRNA))

a2 <- setdiff(de.targets$x2_LIHCGene_TargetScan031000_samples_clustNo_2$miRNA,
              union(de.targets$x2_LIHCGene_TargetScan031000_samples_clustNo_1$miRNA,
                    de.targets$x2_LIHCGene_TargetScan031000_samples_clustNo_3$miRNA))


a3 <- setdiff(de.targets$x2_LIHCGene_TargetScan031000_samples_clustNo_3$miRNA,
              union(de.targets$x2_LIHCGene_TargetScan031000_samples_clustNo_1$miRNA,
                    de.targets$x2_LIHCGene_TargetScan031000_samples_clustNo_2$miRNA))


write.csv(a1, "../test_cases/LIHC/Output/DE_prioritzed_TargetScan_Unique_1.csv")
write.csv(a2, "../test_cases/LIHC/Output/DE_prioritzed_TargetScan_Unique_2.csv")
write.csv(a3, "../test_cases/LIHC/Output/DE_prioritzed_TargetScan_Unique_3.csv")




de.targets <- Reduce(rbind, de.targets)



write.csv(de.targets, "../test_cases/LIHC/Output/DE_prioritzed_TargetScan.csv")







pathways.sets   <- readRDS('../../Data/GeneSets/MSigDB.RDS')
mir.sets        <- readRDS('../../Data/preprocessed/NORMALIZED_MIRSETS_TargetScan03.rds')


miRNA.reporter <- function(miRNA,
                           target_pathway_list,
                           DE_genes_list = NA,
                           miRNAs_list,
                           global_pathways_list){
    
    mir_targets   <- miRNAs_list[[miRNA]]
    path_targets0 <- global_pathways_list[target_pathway_list]
    path_targets  <- Reduce(union,path_targets0)
    path_targets  <- path_targets[path_targets %in% mir_targets]
    
    path_targets2 <- lapply(path_targets0, function(X){
                        X[X %in% path_targets]})
    
    path_targets2 <- path_targets2[lapply(path_targets2, length) > 0]
    path_targets2 <- lapply(path_targets2, function(X){
        path_targets %in% X})
    
    path_targets3 <- Reduce(rbind,path_targets2)
    path_targets3 <- as.data.frame(path_targets3)
    
    rownames(path_targets3) <- names(path_targets2)
    colnames(path_targets3) <- path_targets
    
    sym_targets <- clusterProfiler::bitr(path_targets,fromType = "ENTREZID",
                          toType = "SYMBOL",
                          OrgDb = org.Hs.eg.db)
    
    sym_targets <- sym_targets[!duplicated(sym_targets$ENTREZID),]
    sym_targets <- sym_targets$SYMBOL
    temp.vec    <- apply(path_targets3,1 ,any)
    
    
    path_targets3           <- path_targets3[temp.vec,]
    path_targets4           <- apply(path_targets3, 2,as.numeric)
    colnames(path_targets4) <- sym_targets
    rownames(path_targets4) <- rownames(path_targets3)
    
    if(!is.na(DE_genes_list)){

        mir_tars <- data.frame("SYMBOL" = colnames(path_targets4))
        mir_tars <- left_join(mir_tars,DE_genes_list)
        
        mir_tars[is.na(mir_tars$Direction),]$Direction <- "NONE"
        
        col.labs  <-  levels(factor(mir_tars$Direction))
        mir_tars  <- diag(as.numeric(factor(mir_tars$Direction)))
        mir_tars2 <- path_targets4 %*% mir_tars
        colnames(mir_tars2) <-  colnames(path_targets4)    
    }
    
    return(path_targets4)
    
    
    
    
}

target_pathway_list <- read.csv("~/Desktop/Research/alan_project/Pipeline/temp/Pathways_cluster_louvain.csv",row.names = "X")

target_pathway_list <- target_pathway_list[target_pathway_list$cluster == 3,1]
target_pathway_list <- as.character(target_pathway_list)


library(org.Hs.eg.db)

DE_genes_list       <- readRDS("../Data/LIHC_DE_genes_2021.RDS")
DE_genes_sym        <- clusterProfiler::bitr(rownames(DE_genes_list),
                                             fromType = "ENSEMBL",
                                             toType = "SYMBOL",
                                             OrgDb = org.Hs.eg.db)

DE_genes_list       <- rownames_to_column(DE_genes_list,var = "ENSEMBL")
DE_genes_list       <- left_join(DE_genes_list,DE_genes_sym)

DE_genes_list       <- DE_genes_list %>% dplyr::mutate(Direction = logFC/abs(logFC),
                                     Direction = factor(Direction, 
                                     c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                                     Direction = as.character(Direction))

DE_genes_list$Direction[DE_genes_list$adj.P.Val > 0.05 | abs(DE_genes_list$logFC) < log2(1.5)] = 'NONE'

table(DE_genes_list$Direction)


de.targets
miRNA   <- "hsa-miR-128-3p"	





temp <- miRNA.reporter(miRNA = miRNA,
               target_pathway_list =  target_pathway_list,
               DE_genes_list = DE_genes_list, 
               miRNAs_list = mir.sets,
               global_pathways_list = pathways.sets)

library(ComplexHeatmap)

library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#ffffff", "#999999", "#ef8a62"))
col_fun(seq(-1, 1))


temp2 <- temp[,colSums(temp)>0]
rownames(temp2) <- gsub("Pathway.","",rownames(temp2))

Heatmap(temp2,
        show_column_dend = F,
        show_row_dend = F,
        col = col_fun,
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),
        show_heatmap_legend = F,
        column_title = paste0("Targets of ",miRNA, " in the third cluster (TargetScan)"),
        width = ncol(temp2)*unit(2, "mm"), 
        height = nrow(temp2)*unit(3, "mm"))

sym_targets <- clusterProfiler::bitr(colnames(temp),fromType = "ENTREZID",
                                     toType = "SYMBOL",
                                     OrgDb = org.Hs.eg.db)
