rm(list = ls())
source('01-DifferentialPathwayAnalysis.R')
source('02-MappingPathwaysClusters.R')
source('03-miRNAPathwayEnrichment.R')
source('04-miRNAPrioritization.R')
source('05-miRNAPathwayCorrelation.R')


library(ggplot2)
library(ggfortify)
covariates <- read.csv('../Data/TCGA-LIHC-COV.csv', row.names = 1)
out.dir   = '../test_cases/LIHC/Output/'
data.dir  = '../test_cases/LIHC/Data/'





output0  <- readRDS(paste0(data.dir,"Pathways_object.RDS"))
gene.res <- readRDS(paste0(data.dir,"gene_residuals.RDS"))


PC.path <- prcomp(t(output0$PathwayResiduals), scale.=T, center = T)
PC.gene <- prcomp(t(gene.res), scale.=T, center = T)



# 
# pca <- prcomp(t(TMP.EXP))

p1 <- autoplot(PC.path,
               data = covariates,
               colour = "shortLetterCode",
               size   = 4) 
p1 <- p1 + theme_bw(base_size = 20) + ggtitle("PCA of pathways") + labs(fill = "Tissue type")



p2 <- autoplot(PC.gene, 
               data = covariates,
               colour = "shortLetterCode",
               size   = 4)
p2 <- p2 + theme_bw(base_size = 20) + ggtitle("PCA of genes") + labs(fill = "Tissue type") 

ggsave(paste0(out.dir,"Figures/pca_pathways.pdf"),p1, width = 12, height = 10)
ggsave(paste0(out.dir,"Figures/pca_genes.pdf"),p2, width = 12, height = 10)






path.table         <- output0$DEP
path.table         <- rownames_to_column(path.table,"Pathway")
path.table$Pathway <- gsub("Pathway.","",path.table$Pathway)
path.table$Pathway <- gsub("^([a-zA-Z]*)_","\\1: ",path.table$Pathway)
path.table$Pathway <- gsub("_"," ",path.table$Pathway)

path.table$Direction <- "NONE"
path.table[path.table$logFC<0,]$Direction <- "DOWN"
path.table[path.table$logFC>0,]$Direction <- "UP"

path.table[path.table$adj.P.Val>0.00000001,]$Direction <- "NONE"
path.table <- path.table %>% mutate_if(., is.numeric, signif, digits =3)
write.csv(path.table,paste0(out.dir,"DE_pathways.csv"),row.names = F)



enrichTable <- read.csv(paste0(out.dir,"Enrichment_pathways.csv"))


enrichTable <- as.character(enrichTable[enrichTable$p.adjust< 0.01,]$ID)


intersect(path.table[path.table$Direction != "NONE",]$Pathway,enrichTable)
