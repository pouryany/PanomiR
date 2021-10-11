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
out.dir   = '../test_cases/LIHC/Output/'
data.dir  = '../test_cases/LIHC/Data/'

output0 <- DifferentialPathwayAnalysis(genes.counts,
                                       pathways,
                                       covariates,
                                       condition,
                                       adjust.covars='plate')






test.labels  <-  performDPAwithLabelPertubation(genes.counts,
                                                pathways,
                                                covariates,
                                                condition,
                                                adjust.covars='plate',
                                                samp.rate=1000,
                                                pval.test=F)


test.labels2 <- test.labels$iterations


mean(sapply(test.labels2,function(X){nrow(X[X$adj.P.Val < 0.01,])}))






test.labels3  <-  performDPAwithGenePertubation2(genes.counts,
                                                pathways,
                                                covariates,
                                                condition,
                                                adjust.covars='plate',
                                                samp.rate=1000,
                                                pval.test=F,
                                                numcores = 14)


#saveRDS(test.labels3,"../test_cases/LIHC/Data/Pathway_Randomization_Results.RDS")


all_ps <- sapply(test.labels3,function(X){return(tibble(pval = X$adj.P.Val))})
all_ps <- Reduce(cbind,all_ps)


colnames(all_ps) <- paste0("Permutation_",1:1000)

all_ps <- data.frame(all_ps)


all_ps$average  <- apply(all_ps, 1, mean)
all_ps$observed <- output0$DEP$adj.P.Val

library(reshape2)


all_ps2 <- reshape2::melt(all_ps)
all_ps2$dataset <- "sampled"
all_ps2[all_ps2$variable == "observed",]$dataset <- "observed"
all_ps2[all_ps2$variable == "average",]$dataset <- "average"

all_ps2$dataset <- factor(all_ps2$dataset)



ggplot(all_ps2, aes(x=value, color=dataset,
                    group=variable, alpha = dataset)) +
  stat_ecdf(data = all_ps2, size = .5) +
  theme_bw(base_size = 20) +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values = c("#d8b365", "#e34a33", "#636363" )) +
  scale_alpha_discrete(range = c(0.9, 0.9, 0.2)) + 
  geom_vline(xintercept=0.1,
             colour="grey",
             linetype = "dashed") 




sd(sapply(test.labels3,function(X){nrow(X[X$adj.P.Val < 0.01,])}))

nrow(output0$DEP[output0$DEP$adj.P.Val < 0.01,])
pnorm(825,693.785,32.90945,lower.tail = F)
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
                                            out.dir="../temp22",
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
mir.selection   <- names(mir.sets)



for (mirs in mir.sets.list){
    
    tag       <- tail(unlist(stringr::str_split(mirs,pattern = "_")),1)
    mir.sets  <- readRDS(mirs)
    name.tag  <- paste0("LIHCGenesLIHCMirsENRICHMENT_",tag)
    
    print(paste0("performing: ", tag))
    enriches0 <- miRNAPathwayEnrichment(mir.sets,
                                        pathways.sets,
                                        genes.selection = genes.selection,
                                        mir.selection = mir.selection,
                                        save.RDS.name = name.tag,
                                        out.dir= data.dir)
    
}



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
        output2 <- PrioritizeMicroRNA(enriches0,
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
        output2 <- PrioritizeMicroRNA(enriches0,
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

