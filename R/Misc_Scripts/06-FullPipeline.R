rm(list = ls())
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


if(!dir.exists(out.dir))
    dir.create(out.dir,recursive = T)
if(!dir.exists(data.dir))
    dir.create(data.dir,recursive = T)

path.table   <- output0$DEP
path.table   <- rownames_to_column(path.table,"Pathway")
path.table$Pathway <- gsub("Pathway.","",path.table$Pathway)
path.table$Pathway <- gsub("^([a-zA-Z]*)_","\\1: ",path.table$Pathway)
path.table$Pathway <- gsub("_"," ",path.table$Pathway)

path.table$Direction <- "NONE"
path.table[path.table$logFC<0,]$Direction <- "DOWN"
path.table[path.table$logFC>0,]$Direction <- "UP"

path.table[path.table$adj.P.Val>0.001,]$Direction <- "NONE"
path.table <- path.table %>% mutate_if(., is.numeric, signif, digits =3)
write.csv(path.table,paste0(out.dir,"DE_pathways.csv"),row.names = F)





de.paths <- output0$DEP
saveRDS(de.paths,paste0(out.dir,"LIHC_DEP.RDS"))
de.paths <- readRDS(paste0(out.dir,"LIHC_DEP.RDS"))

# pathway.summaries <- output0$pathwaySummaryStats
# pathway.summaries <- readRDS('../Data/TCGA-LIHC-x2-all-PathwaySummaryStats-NEW.RDS')

pcxn <- readRDS('../Data/GeneSets/improved_PCxN_MSigDB.RDS')

pathway.clusters <- MappingPathwaysClusters(pcxn = pcxn, 
                                            de.paths = de.paths[1:300,],
                                            out.dir=out.dir,
                                            subplot = F,
                                            prefix='LIHC_x2_all_top200_',
                                            cor.thresh = 0.1)



mir.sets        <- readRDS('../Data/preprocessed/NORMALIZED_MIRSETS.rds')
pathways.sets   <- readRDS('../Data/GeneSets/MSigDB.RDS')
genes.selection <- rownames(genes.counts)
mirna.counts    <- readRDS('../Data/TCGA-LIHC-DEmiRNA.RDS')
mir.selection   <- rownames(mirna.counts)


# Set up default enrichment results if detail are not provided. 
# Set up a selection of enrichment dataset. 


enriches0 <- miRNAPathwayEnrichment(mir.sets,
                                    pathways.sets,
                                    genes.selection = genes.selection,
                                    mir.selection = mir.selection,
                                    save.RDS.name = 'LIHCGenesLIHCMirsENRICHMENT.RDS',
                                    out.dir= data.dir)


top.clusters <- names(sort(table(pathway.clusters$cluster),decreasing = T))[1:2]
top.clusters <- pathway.clusters[pathway.clusters$cluster %in% top.clusters,]

method <- c('AggInv')
output <- miRNAPrioritization(enriches0,
                              top.clusters,
                              method,
        ß                      out.dir=paste0(out.dir, 'Prioritization/'),
                              data.dir=data.dir,
                              samp.rate=1000,
                              prefix='x2_all_LIHCGenesLIHCMirs_NEW_',
                              save.jack.knife=T,
                              save.csv=T)
temp <- output$Cluster1
# output1 <- miRNAPathwayCorrelation(pathway.summaries,
#                                    pathway.clusters,
#                                    mirna.counts,
#                                    covariates,
#                                    condition,
#                                    prime.cond = 'TP',
#                                    compare.cond = 'NT',
#                                    adjust.covars = 'plate')
