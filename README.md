PanomiR Walkthrough
-------------------

PanomiR is a package to detect miRNAs that target groups of pathways
from gene expression data. Here, we describe the work flow of PanomiR
along with different steps of the package. The working directory should
be set to the `R` folder.

The following are initial set for the program. The main input at this
stage are as follows:

1.  A table of pathway membership. In each row there is a name of a
    pathway and a gene that belongs to it.
2.  A gene expression dataset.
3.  A table of covariates.
4.  A list of covariates that you would like to adjust for.
5.  Output and datadirectories.

Here we use an example from TCGA Liver Hepatocellular Carcinoma (LIHC)
to demonstrate the utility of PanomiR.

``` r
source('01-DifferentialPathwayAnalysis.R')
source('02-MappingPathwaysClusters.R')
source('03-miRNAPathwayEnrichment.R')
source('04-miRNAPrioritization.R')
source('05-miRNAPathwayCorrelation.R')


pathways     <- readRDS('../Data/preprocessed/MSigDBPathGeneTab.RDS')
genes.counts <- readRDS('../Data/LIHC_gene.RDS')
covariates   <- read.csv('../Data/TCGA-LIHC-COV.csv', row.names = 1)



condition = 'shortLetterCode'
out.dir   = '../test_cases2/LIHC/Output/'
data.dir  = '../test_cases2/LIHC/Data/'

output0 <- DifferentialPathwayAnalysis(genes.counts,
                                       pathways,
                                       covariates,
                                       condition,
                                       adjust.covars='plate')
```

Next step is to Map differentially expressed pathways to a precalculated
network of pathway correlations to find modules of DE-pathways. The
function `MappingPathwayClusters` can take different graph clustering
algorithms. Here we use a for loop to determine readouts for different
clustering methods.

``` r
de.paths <- output0$DEP

# Precalculated
pcxn      <- readRDS('../Data/GeneSets/improved_PCxN_MSigDB.RDS')
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
```

Next for each miRNA background and for each clustering algorithm, the
miRNA targeting scores are calculated. Choice of appropriate algorithms
are left to users.

``` r
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
```

For miRNA prioritization it is important to have a precalculated list of
miRNA-pathway association for all possible pairs. While in the
distributed package we will provide some tables. Here, I will
demonstrate how to generate such tables using miRNA-target lists.

``` r
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
```

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
