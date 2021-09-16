rm(list = ls())
source('01-DifferentialPathwayAnalysis.R')
source('02-MappingPathwaysClusters.R')
source('03-miRNAPathwayEnrichment.R')
source('04-miRNAPrioritization.R')
source('05-miRNAPathwayCorrelation.R')


library(ggplot2)

pathways     <- readRDS('../Data/preprocessed/MSigDBPathGeneTab.RDS')
genes.counts <- readRDS('../Data/LIHC_gene.RDS')
covariates <- read.csv('../Data/TCGA-LIHC-COV.csv', row.names = 1)
condition = 'shortLetterCode'
out.dir   = '../test_cases2/LIHC/Output/'
data.dir  = '../test_cases2/LIHC/Data/'

output0 <- DifferentialPathwayAnalysis(genes.counts,
                                       pathways,
                                       covariates,
                                       condition,
                                       adjust.covars='plate')


if(!dir.exists(out.dir))
    dir.create(out.dir,recursive = T)
if(!dir.exists(data.dir))
    dir.create(data.dir,recursive = T)


saveRDS(output0,paste0(data.dir,"Pathways_object.RDS"))





de.paths <- output0$DEP
saveRDS(de.paths,paste0(out.dir,"LIHC_DEP.RDS"))
de.paths0 <- readRDS(paste0(data.dir,"Pathways_object.RDS"))
de.paths  <- de.paths0$DEP
# pathway.summaries <- output0$pathwaySummaryStats
# pathway.summaries <- readRDS('../Data/TCGA-LIHC-x2-all-PathwaySummaryStats-NEW.RDS')

pcxn <- readRDS('../Data/GeneSets/improved_PCxN_MSigDB.RDS')

pathway.clusters <- MappingPathwaysClusters(pcxn = pcxn, 
                                            de.paths = de.paths[1:300,],
                                            out.dir=out.dir,
                                            subplot = F,
                                            prefix='LIHC_x2_all_top200_',
                                            cor.thresh = 0.1)


func_list <- c("cluster_edge_betweenness",
               "cluster_infomap",
               "cluster_fast_greedy",
               "cluster_louvain")

pathway.clusters <- list()
for(func in func_list){
    
    temp.clusters <- MappingPathwaysClusters(pcxn = pcxn, 
                                            de.paths = de.paths[1:300,],
                                            out.dir= out.dir,
                                            subplot = F, 
                                            top.paths = 200,
                                            prefix= paste0('top200_',func),
                                            cor.thresh = 0.1,
                                            clust.fn = get(func),
                                            save.csv.name = paste0("Pathways_",
                                                                    func,
                                                                    ".csv"))
    temp.clusters$method     <- func
    pathway.clusters[[func]] <- temp.clusters
    
}



Reduce(rbind,pathway.clusters)


mir.sets        <- readRDS('../Data/preprocessed/NORMALIZED_MIRSETS.rds')
#mir.sets        <- readRDS('../../Data/preprocessed/NORMALIZED_MIRSETS_TargetScan01.rds')
mir.sets.list    <- list.files("../Data/preprocessed/",
                              pattern = "TargetScan",
                              full.names = T)


pathways.sets   <- readRDS('../Data/GeneSets/MSigDB.RDS')
genes.selection <- rownames(genes.counts)
mirna.counts    <- readRDS('../Data/TCGA-LIHC-miRNAs_residuals.RDS')
# mir.selection   <- rownames(mirna.counts)
mir.selection   <- names(mir.sets)

# Set up default enrichment results if detail are not provided. 


# 
# 
# 
# enriches0 <- miRNAPathwayEnrichment(mir.sets,
#                                     pathways.sets,
#                                     genes.selection = genes.selection,
#                                     mir.selection = mir.selection,
#                                     save.RDS.name = 'LIHCGenesLIHCMirsENRICHMENT_Tarbase.RDS',
#                                     out.dir= data.dir)
# 





# Set up a selection of enrichment dataset. 


for (mirs in mir.sets.list){
    
  
    tag       <- tail(unlist(stringr::str_split(mirs,pattern = "_")),1)
    mir.sets  <- readRDS(mirs)
    name.tag  <- paste0("LIHCGenesLIHCMirsENRICHMENT_",tag)
    
    mir.selection2 <- names(mir.sets)
    
    print(paste0("performing: ", tag))
    enriches0 <- miRNAPathwayEnrichment(mir.sets,
                                        pathways.sets,
                                        genes.selection = genes.selection,
                                        mir.selection = mir.selection2,
                                        save.RDS.name = name.tag,
                                        out.dir= data.dir)
    
}




# 
#  enriches0 <- miRNAPathwayEnrichment(mir.sets,
#                                      pathways.sets,
#                                      genes.selection = genes.selection,
#                                      mir.selection = mir.selection,
#                                      save.RDS.name = 'LIHCGenesLIHCMirsENRICHMENT_Tarbase.RDS',
#                                      out.dir= data.dir)
#  
#  
#  

  enriches0 <- readRDS(paste0(data.dir,"LIHCGenesLIHCMirsENRICHMENT_Tarbase.RDS"))
# enriches0 <- readRDS(paste0(data.dir,"LIHCGenesLIHCMirsENRICHMENT_TargetScan01.RDS"))

# top.clusters <- names(sort(table(pathway.clusters$cluster),decreasing = T))[1:2]

# top.clusters <- names(sort(table(pathway.clusters$cluster_louvain$cluster),decreasing = T))[1:2]
# top.clusters <- pathway.clusters$cluster_louvain[pathway.clusters$cluster_louvain$cluster %in% top.clusters,]
# 
# # method <- c('AggInv')
# # output <- miRNAPrioritization(enriches0,
# #                               top.clusters,
# #                               method,
# #                               out.dir=paste0(out.dir, 'Prioritization/'),
# #                               data.dir=data.dir,
# #                               samp.rate=1000,
# #                               prefix='x2_all_LIHCGenesLIHCMirs_NEW_',
# #                               save.jack.knife=T,
# #                               save.csv=T)
# # 
# # 
# 
# method <- c('AggInv')
# output2 <- miRNAPrioritization2(enriches0,
#                               top.clusters,
#                               method,
#                               out.dir=paste0(out.dir, 'Prioritization/'),
#                               data.dir=data.dir,
#                               samp.rate=2000,
#                               prefix='x2_all_LIHCGenesLIHCMirs_NEW_',
#                               save.jack.knife=T,
#                               save.csv=T,num.cores = 8)
# 


  enriches0 <- readRDS(paste0(data.dir,"LIHCGenesLIHCMirsENRICHMENT_Tarbase.RDS"))
  
func_list <- c("cluster_edge_betweenness",
               "cluster_infomap",
               "cluster_fast_greedy",
               "cluster_louvain")



for(func in func_list){
    
    method <- c('AggInv')
    
    top.clusters <- pathway.clusters[[func]]
    enriches0    <- readRDS(paste0(data.dir,"LIHCGenesLIHCMirsENRICHMENT_Tarbase.RDS"))
    
        print(paste0("performing: ", func))
        output2 <- miRNAPrioritization2(enriches0,
                                        top.clusters,
                                        method,
                                        out.dir=paste0(out.dir,func,
                                                       '_Prioritization_',
                                                       "Tarbase",
                                                       '/'),
                                        data.dir=data.dir,
                                        samp.rate=1000,
                                        prefix=paste0('x2_LIHCGene_',"Tarbase"),
                                        save.jack.knife=F,
                                        save.csv=T,
                                        num.cores = 8,
                                        top.clust=3)

    
    
}






for(func in func_list){
    
    method       <- c('AggInv')
    top.clusters <- pathway.clusters[[func]]
    
    for (mirs in mir.sets.list){
        
        tag       <- tail(unlist(stringr::str_split(mirs,pattern = "_")),1)
        mir.sets  <- readRDS(mirs)
        name.tag  <- paste0("LIHCGenesLIHCMirsENRICHMENT_",tag)
        
        tag       <- gsub(".rds","",tag)
        
        enriches0 <- readRDS(paste0(data.dir,name.tag))
        
        print(paste0("performing: ", tag))
        output2 <- miRNAPrioritization2(enriches0,
                                        top.clusters,
                                        method,
                                        out.dir=paste0(out.dir,func,
                                                       '_Prioritization_',
                                                       tag,
                                                       '/'),
                                        data.dir=data.dir,
                                        samp.rate=1000,
                                        prefix=paste0('x2_LIHCGene_',tag),
                                        save.jack.knife=F,
                                        save.csv=T,
                                        num.cores = 8,
                                        top.clust=3)
        
    }

    
}




# Testing the effect of PCxN mapping on clustering 




pathway.clusters <- list()
func <- "cluster_louvain"
for(n_paths in c(150,200,250,300,350,400,450)){
  
  temp.clusters <- MappingPathwaysClusters(pcxn = pcxn, 
                                           de.paths = de.paths[1:500,],
                                           out.dir="../temp22",
                                           subplot = F, 
                                           top.paths = n_paths,
                                           prefix= "",
                                           cor.thresh = 0.1,
                                           clust.fn = get(func),
                                           save.csv.name = NULL)
  temp.clusters$method        <- func
  temp.clusters$n_paths       <- n_paths
  pathway.clusters[[n_paths]] <- temp.clusters
  print(paste0("Done ", n_paths))
  
}



path.reps <- Reduce(rbind,pathway.clusters)

path.reps$tag <- paste0(path.reps$cluster,"_",path.reps$n_paths)
path.reps     <- path.reps[path.reps$cluster %in% 1:3,]

tagz <- unique(path.reps$tag)

tagz.mat <- matrix(0, length(tagz), length(tagz))
rownames(tagz.mat) <- colnames(tagz.mat) <- tagz

jaccard.ind <- function(set1,set2){
  uns <- length(union(set1,set2))
  ovr <- length(intersect(set1,set2))
  jac <- ovr/uns
}



for(set1 in tagz){
  for(set2 in tagz){
    set1.paths <- path.reps[path.reps$tag == set1,]$Pathway
    set2.paths <- path.reps[path.reps$tag == set2,]$Pathway
    
    tagz.mat[set1,set2] <- jaccard.ind(set1.paths,set2.paths)
    
  }
}

heatmap(tagz.mat)

pheatmap::pheatmap(tagz.mat,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   color=colorRampPalette(c("#f0f0f0","#de2d26"))(50))

?pheatmap
# 
# 
# 
# temp2<- output2$Cluster1
# 
# temp <- output$Cluster1
# 
# 
# mirna.counts2 <- mirna.counts
# colnames(mirna.counts2) <- rownames(covariates)
# 
# 
# 
# mirna.data <- read.csv('../Data/LIHC_miRseq_mature_RPM.txt', sep='\t')
# rownames(mirna.data)<- sapply(as.character(mirna.data$Gene), function(X){
#     return(substr(X, 1, nchar(X)-13))
# })
# mirna.data <- mirna.data %>% dplyr::select(., -c('Gene'))
# colnames(mirna.data) <- gsub('\\.', '-', colnames(mirna.data))
# mirna.data[is.na(mirna.data)] <- 0
# 
# 
# colnames(mirna.data) <- gsub("\\.","-",colnames(mirna.data))
# mirna.data           <- mirna.data[rownames(mirna.counts),colnames(mirna.counts)]
# colnames(mirna.data) <- rownames(covariates)
# mirna.data           <- as.matrix(mirna.data)
# 
# output1 <- miRNAPathwayCorrelation(de.paths0$PathwayResiduals,
#                                     top.clusters[,c(1,2)],
#                                     mirna.counts2,
#                                     covariates,
#                                     condition,
#                                     prime.cond = 'TP',
#                                     compare.cond = 'NT',
#                                     top.clust = 2,
#                                     num.cores = 8, cor.method = "pearson")
# 
# 
# temp2 <- output1$Cluster1
# 
# 
# 
# 
# tarb1 <- tarb[order(tarb$AggInv_pval),]
# targ1 <- targ[order(targ$AggInv_pval),]
# 
# tarb1 <- tarb1[1:100,]$x
# targ1 <- targ1[1:100,]$x
# 
# intersect(tarb1,targ1)
# 
# DE_genes <- rownames(`TCGA-LIHC-DEmiRNAs`[1:100,])
# 
# 
# DE_genes[DE_genes %in% targ$x]
# targ2 <- targ[grep(paste(DE_genes,collapse = "|"), targ$x),]
# 
# intersect(DE_genes,targ1)
