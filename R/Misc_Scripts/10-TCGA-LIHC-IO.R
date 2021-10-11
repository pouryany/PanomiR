library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
source('revised_functions.R')

project <- 'LIHC'

query <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

query2 <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")

# GDCdownload(query,directory = '../GDCdata')
# GDCdownload(query2, directory = '../GDCdata')

htseq.summarized.exp <- GDCprepare(query, directory = '../GDCdata')

htseq.samplesDown <- getResults(query,cols = c("cases"))

htseq.SmTP <- TCGAquery_SampleTypes(barcode = htseq.samplesDown,
                                  typesample = "TP")
htseq.SmNT <- TCGAquery_SampleTypes(barcode = htseq.samplesDown,
                                  typesample = "NT")
# 
# mirna.samplesDown <- getResults(query2,cols = c("cases"))
# 
# mirna.SmTP <- TCGAquery_SampleTypes(barcode = mirna.samplesDown,
#                                     typesample = "TP")
# mirna.SmNT <- TCGAquery_SampleTypes(barcode = mirna.samplesDown,
#                                     typesample = "NT")
# 
# htseq.TP.samples <- sort(htseq.SmTP[substr(htseq.SmTP, 1, 12) %in%  substr(mirna.SmTP, 1, 12)])
# htseq.NT.samples <- sort(htseq.SmNT[substr(htseq.SmNT, 1, 12) %in%  substr(mirna.SmNT, 1, 12)])
# 
# mirna.TP.samples <- sort(mirna.SmTP[substr(mirna.SmTP, 1, 12) %in%  substr(htseq.SmTP, 1, 12)])
# mirna.NT.samples <- sort(mirna.SmNT[substr(mirna.SmNT, 1, 12) %in%  substr(htseq.SmNT, 1, 12)])
# 

dataPrep <- TCGAanalyze_Preprocessing(object = htseq.summarized.exp, 
                          cor.cut = 0.6,
                          datatype = "HTSeq - Counts")

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

dataNorm2 <- TCGAanalyze_Normalization(tabDF = dataNorm,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                   method = "quantile", 
                                   qnt.cut =  0.25)

filter <- (rowSums(dataFilt > 1) >= floor(nrow(covariates)/2))
sum(filter)

dataFilt <- dataFilt[filter,rownames(covariates)]

metadata <- colData(htseq.summarized.exp)
metadata <- metadata[colnames(dataFilt),]
metadata <- as.data.frame(metadata)
metadata.treament <- dplyr::select(metadata, c('treatments'))
metadata <- dplyr::select(metadata, -c('treatments'))

extension <- get_IDs(dataFilt)
extension <- dplyr::select(extension, -c(patient, sample))
metadata <- merge(metadata, extension, by='barcode')
rownames(metadata) <- metadata$barcode
metadata <- metadata[colnames(dataFilt),]
metadata$mirna_sample <- c(mirna.TP.samples, mirna.NT.samples)

# Error in plot.window(...) : need finite 'ylim' values; Error in plotting

# batch.corrected <- TCGAbatch_Correction(tabDF = dataFilt, batch.factor = 'Plate')

# Added plotting options; plot=FALSE to skip plotting
batch.corrected <- TCGAbatch_Correction_revised(tabDF = dataFilt, batch.factor = 'Plate', plot = FALSE)

IDs <- get_IDs(dataFilt)
batch.corrected.limma <- removeBatchEffect(dataFilt, IDs$plate)

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,htseq.TP.samples],
                            mat2 = dataFilt[,htseq.NT.samples],
                            batch.factors = 'Plate',
                            pipeline="limma",
                            Cond1type = "Tumor",
                            Cond2type = "Normal",
                            fdr.cut = 1,
                            log.trans = T,
                            logFC.cut = 0,
                            voom = T,
                            ClinicalDF = data.frame())


dataDEGs.batch <- TCGAanalyze_DEA(mat1 = batch.corrected[,dataSmTP],
                            mat2 = batch.corrected[,dataSmNT],
                            pipeline="limma",
                            Cond1type = "Tumor",
                            Cond2type = "Normal",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            ClinicalDF = data.frame())

dataDEGs.batch.limma <- TCGAanalyze_DEA(mat1 = batch.corrected.limma[,dataSmTP],
                                  mat2 = batch.corrected.limma[,dataSmNT],
                                  pipeline="limma",
                                  Cond1type = "Tumor",
                                  Cond2type = "Normal",
                                  fdr.cut = 0.01 ,
                                  logFC.cut = 1,
                                  ClinicalDF = data.frame())

write.csv(metadata, '../Data/covariates.csv')

saveRDS(dataFilt, file='../Data/TCGA-LIHC-GENE-COUNTS.RDS')
saveRDS(rownames(dataFilt), file='../Data/LIHC-GENES.RDS')


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

