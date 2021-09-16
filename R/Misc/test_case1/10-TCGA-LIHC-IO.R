rm(list = ls())
library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
source('00-revised_functions.R')

project   = 'LIHC'
data.dir  = '../test_cases/LIHC/Data/'
out.dir   = '../test_cases/LIHC/Output/'



if(FALSE){
    
    query <- GDCquery(project = "TCGA-LIHC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")
    # 
    # query2 <- GDCquery(project = "TCGA-LIHC",
    #                   data.category = "Transcriptome Profiling",
    #                   data.type = "miRNA Expression Quantification")
    
    # GDCdownload(query,directory = '../GDCdata')
    # GDCdownload(query2, directory = '../GDCdata')
    
    htseq.summarized.exp <- GDCprepare(query, directory = '../GDCdata')
    
    htseq.samplesDown <- getResults(query,cols = c("cases"))
    # 
    htseq.SmTP <- TCGAquery_SampleTypes(barcode = htseq.samplesDown,
                                        typesample = "TP")
    htseq.SmNT <- TCGAquery_SampleTypes(barcode = htseq.samplesDown,
                                        typesample = "NT")
    # # 
    mirna.samplesDown <- getResults(query2,cols = c("cases"))
    # 
    mirna.SmTP <- TCGAquery_SampleTypes(barcode = mirna.samplesDown,
                                        typesample = "TP")
    mirna.SmNT <- TCGAquery_SampleTypes(barcode = mirna.samplesDown,
                                        typesample = "NT")
    # 
    htseq.TP.samples <- sort(htseq.SmTP[substr(htseq.SmTP, 1, 12) %in%  substr(mirna.SmTP, 1, 12)])
    htseq.NT.samples <- sort(htseq.SmNT[substr(htseq.SmNT, 1, 12) %in%  substr(mirna.SmNT, 1, 12)])
    # 
    mirna.TP.samples <- sort(mirna.SmTP[substr(mirna.SmTP, 1, 12) %in%  substr(htseq.SmTP, 1, 12)])
    mirna.NT.samples <- sort(mirna.SmNT[substr(mirna.SmNT, 1, 12) %in%  substr(htseq.SmNT, 1, 12)])
    # 
    
    # covariates of samples with concordant miRNA/gene expression
    covariates <- c(htseq.NT.samples,htseq.TP.samples)
    
    
    saveRDS(htseq.summarized.exp,paste0(data.dir,"GeneSummarizedExperiment.RDS"))
    
}


htseq.summarized.exp <- readRDS(paste0(data.dir,"GeneSummarizedExperiment.RDS"))


dataPrep <- TCGAanalyze_Preprocessing(object = htseq.summarized.exp, 
                          cor.cut = 0.6,
                          datatype = "HTSeq - Counts")

# dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                       geneInfo = geneInfoHT,
#                                       method = "gcContent")
# 
# dataNorm2 <- TCGAanalyze_Normalization(tabDF = dataNorm,
#                                       geneInfo = geneInfoHT,
#                                       method = "geneLength")
# 
# dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                    method = "quantile", 
#                                    qnt.cut =  0.25)
# 
# 
# filter <- (rowSums(dataFilt > 1) >= floor(ncol(dataNorm)/5))
# sum(filter)
# 
# dataFilt <- dataFilt[filter,(covariates)]
# 
# metadata <- colData(htseq.summarized.exp)
# metadata <- metadata[colnames(dataFilt),]
# metadata <- as.data.frame(metadata)
# metadata.treament <- dplyr::select(metadata, c('treatments'))
# metadata <- dplyr::select(metadata, -c('treatments'))
# 
# extension <- get_IDs(dataFilt)
# extension <- dplyr::select(extension, -c(patient, sample))
# metadata <- merge(metadata, extension, by='barcode')
# rownames(metadata) <- metadata$barcode
# metadata <- metadata[colnames(dataFilt),]
# metadata$mirna_sample <- c(mirna.TP.samples, mirna.NT.samples)
# # bad.batch <- names(which(table(metadata$plate) == 1))
# # metadata <- metadata[!(metadata$plate %in% bad.batch),]
# # 
# 
# 
# dataFilt <- dataFilt[,rownames(metadata)]
# 
# htseq.NT.samples  <- rownames(metadata[metadata$shortLetterCode == "NT",])
# htseq.TP.samples  <- rownames(metadata[metadata$shortLetterCode == "TP",])
# # Error in plot.window(...) : need finite 'ylim' values; Error in plotting
# 
# # batch.corrected <- TCGAbatch_Correction(tabDF = dataFilt, batch.factor = 'Plate')
# 
# # Added plotting options; plot=FALSE to skip plotting
# batch.corrected <- TCGAbatch_Correction_revised(tabDF = dataFilt, batch.factor = 'Plate', plot = FALSE)
# 
# IDs <- get_IDs(dataFilt)
# batch.corrected.limma <- removeBatchEffect(dataFilt, IDs$plate)
# 
# dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,htseq.TP.samples],
#                             mat2 = dataFilt[,htseq.NT.samples],
#                             batch.factors = 'Plate',
#                             pipeline="limma",
#                             Cond1type = "Tumor",
#                             Cond2type = "Normal",
#                             fdr.cut = 1,
#                             log.trans = T,
#                             logFC.cut = 0,
#                             voom = T,
#                             ClinicalDF = data.frame())
# 
# 
# dataDEGs.batch <- TCGAanalyze_DEA(mat1 = batch.corrected[,htseq.TP.samples],
#                             mat2 = batch.corrected[,htseq.TP.samples],
#                             pipeline="limma",
#                             Cond1type = "Tumor",
#                             Cond2type = "Normal",
#                             fdr.cut = 0.01 ,
#                             logFC.cut = 1,
#                             ClinicalDF = data.frame())
# 
# dataDEGs.batch.limma <- TCGAanalyze_DEA(mat1 = batch.corrected.limma[,htseq.TP.samples],
#                                   mat2 = batch.corrected.limma[,htseq.TP.samples],
#                                   pipeline="limma",
#                                   Cond1type = "Tumor",
#                                   Cond2type = "Normal",
#                                   fdr.cut = 0.01 ,
#                                   logFC.cut = .1,
#                                   ClinicalDF = data.frame())
# 
# write.csv(metadata, '../Data/covariates.csv')
# 
# saveRDS(dataFilt, file='../Data/TCGA-LIHC-GENE-COUNTS.RDS')
# saveRDS(rownames(dataFilt), file='../Data/LIHC-GENES.RDS')
# 
# 
# 
# 

covariates <- read.csv('../Data/TCGA-LIHC-COV.csv', row.names=1)
dataPrep <- dataPrep[,rownames(covariates)]
genes.filter <- rownames(geneInfoHT)
dataPrep <- dataPrep[rownames(dataPrep) %in% genes.filter,]

count.filter <- rowSums(dataPrep < 1) < (ncol(covariates)/2)
dataPrep <- dataPrep[count.filter,]

dataPrep <- getGeneFilteredGeneExprMatrix(dataPrep, 
                                          MIN_GENE_CPM=1, 
                                          MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.2,
                                          verbose = T)
dataPrep <- dataPrep$filteredExprMatrix

genes.data <- geneInfoHT[rownames(dataPrep),]

library(cqn)
dataNorm <- cqn(dataPrep$counts, 
                x=genes.data$gcContent,
                lengths=genes.data$geneLength,
                lengthMethod = 'smooth',
                verbose = F)
dataNorm$E <- dataNorm$y + dataNorm$offset
new.counts <- (2^dataNorm$E)
saveRDS(new.counts,"../test_datasets/LIHC_gene.RDS")
new.counts <- readRDS("../test_datasets/LIHC_gene.RDS")

#Saving Residuals


design <- model.matrix(~0 + plate, data = covariates)

VOOM.WEIGHTS  <- voom(new.counts, design=design, plot=T)
FIT           <- lmFit(log(new.counts), design,
                       weights = VOOM.WEIGHTS$weights)


gene.residuals <- residuals.MArrayLM(FIT,log(new.counts))

saveRDS(gene.residuals,paste0(data.dir,"gene_residuals.RDS"))




design <- model.matrix(~0 + shortLetterCode + plate, data = covariates)

VOOM.WEIGHTS  <- voom(new.counts, design=design, plot=T)
FIT           <- lmFit(log(new.counts), design,
                       weights = VOOM.WEIGHTS$weights)


contrast      <- makeContrasts("shortLetterCodeTP-shortLetterCodeNT",
                               levels=colnames(FIT$coefficients))

FIT.CONTR     <- contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR     <- eBayes(FIT.CONTR)
tT            <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)



saveRDS(tT,"../Data/LIHC_DE_genes_2021.RDS")
tT <- readRDS("../Data/LIHC_DE_genes_2021.RDS")


library(clusterProfiler)
library(org.Hs.eg.db)

GENE.PARAM   <-   clusterProfiler::bitr(genes.filter,
                                  fromType = "ENSEMBL",
                                  toType = "SYMBOL",
                                  OrgDb = org.Hs.eg.db)


DE = lapply(1:dim(contrast)[2], function(i, FIT){
    topTable(FIT, coef=i, number = 50000, confint = T) %>%
        rownameToFirstColumn('ENSEMBL')
}, FIT.CONTR)
names(DE) = c('Tumor - Normal')

DE = DE %>% 
    rbindlist(idcol = 'Comparison') %>%
    left_join(GENE.PARAM %>%
                  dplyr::select(ENSEMBL, SYMBOL) %>%
                  unique()) %>%
    dplyr::mutate(Direction = logFC/abs(logFC),
                  Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                  Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.05 | abs(DE$logFC) < 1] = 'NONE'



saveRDS(DE,"../test_cases/LIHC/Output/DE_genes.RDS")

pathways     <- readRDS('../../Data/preprocessed/MSigDBPathGeneTab.RDS')
pathways     <- pathways[,c(1,3)]
names(pathways) <- c("term","gene")


enrichment.res <- enricher(DE[DE$Direction != "NONE",]$ENSEMBL, pAdjustMethod = "fdr",
         universe = DE$ENSEMBL,minGSSize = 10, maxGSSize =2000 ,TERM2GENE = pathways)


enrichment.res <- enrichment.res@result

saveRDS(enrichment.res,"../test_cases/LIHC/Output/enrichment.RDS")


path.table         <- enrichment.res
path.table$ID <- gsub("Pathway.","",path.table$ID)
path.table$ID <- gsub("^([a-zA-Z]*)_","\\1: ",path.table$ID)
path.table$ID <- gsub("_"," ",path.table$ID)


path.table <- path.table %>% mutate_if(., is.numeric, signif, digits =3)
write.csv(path.table,paste0(out.dir,"Enrichment_pathways.csv"),row.names = F)

nrow(path.table[path.table$p.adjust < 0.01,])


# 
# 
# 
# 
# DE2 = `TCGA-LIHC-DEG` %>%
#         rownameToFirstColumn('ENSEMBL')
# 
# DE2 <- list(DE2)
# names(DE2) = c('Tumor - Normal')
# 
# DE2 = DE2 %>% 
#     rbindlist(idcol = 'Comparison') %>%
#     left_join(GENE.PARAM %>%
#                   dplyr::select(ENSEMBL, SYMBOL) %>%
#                   unique()) %>%
#     dplyr::mutate(Direction = logFC/abs(logFC),
#                   Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
#                   Direction = as.character(Direction))
# DE2$Direction[DE2$adj.P.Val > 0.05 | abs(DE2$logFC) < 1] = 'NONE'
# 
# 


