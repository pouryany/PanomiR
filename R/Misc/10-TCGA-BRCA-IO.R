library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
source('revised_functions.R')

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

query2 <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification")

# GDCdownload(query, directory = '../GDCdata')
# GDCdownload(query2, directory = '../GDCdata')

htseq.summarized.exp <- GDCprepare(query, directory = '../GDCdata')

htseq.samplesDown <- getResults(query,cols = c("cases"))

htseq.SmTP <- TCGAquery_SampleTypes(barcode = htseq.samplesDown,
                                  typesample = "TP")
htseq.SmNT <- TCGAquery_SampleTypes(barcode = htseq.samplesDown,
                                  typesample = "NT")

mirna.samplesDown <- getResults(query2,cols = c("cases"))

mirna.SmTP <- TCGAquery_SampleTypes(barcode = mirna.samplesDown,
                                    typesample = "TP")
mirna.SmNT <- TCGAquery_SampleTypes(barcode = mirna.samplesDown,
                                    typesample = "NT")

mirna.data <- read.csv('../Data/BRCA_miRseq_mature_RPM.txt', sep='\t')
rownames(mirna.data)<- sapply(as.character(mirna.data$Gene), function(X){
  return(substr(X, 1, nchar(X)-13))
})
mirna.data <- mirna.data %>% dplyr::select(., -c('Gene'))
colnames(mirna.data) <- gsub('\\.', '-', colnames(mirna.data))
mirna.data[is.na(mirna.data)] <- 0

mirna.SmTP <- colnames(mirna.data)[colnames(mirna.data) %in% mirna.SmTP] 
mirna.SmNT <- colnames(mirna.data)[colnames(mirna.data) %in% mirna.SmNT] 

htseq.TP.samples <- sort(htseq.SmTP[substr(htseq.SmTP, 1, 12) %in%  substr(mirna.SmTP, 1, 12)])
htseq.NT.samples <- sort(htseq.SmNT[substr(htseq.SmNT, 1, 12) %in%  substr(mirna.SmNT, 1, 12)])

mirna.TP.samples <- sort(mirna.SmTP[substr(mirna.SmTP, 1, 12) %in%  substr(htseq.SmTP, 1, 12)])
mirna.NT.samples <- sort(mirna.SmNT[substr(mirna.SmNT, 1, 12) %in%  substr(htseq.SmNT, 1, 12)])

mirna.filter <- rowSums(mirna.data > 1) >= floor(length(c(mirna.TP.samples, mirna.NT.samples))/2)
mirna.data <- mirna.data[mirna.filter, c(mirna.TP.samples, mirna.NT.samples)]
saveRDS(mirna.data, '../Data/TCGA-BRCA-miRNA-COUNTS.RDS')

dataPrep <- TCGAanalyze_Preprocessing(object = htseq.summarized.exp, 
                          cor.cut = 0.6,
                          datatype = "HTSeq - Counts")

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

dataNorm <- TCGAanalyze_Normalization(tabDF = dataNorm,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                   method = "quantile", 
                                   qnt.cut =  0.25)

filter <- (rowSums(dataFilt > 10) >= floor(length(c(htseq.TP.samples, htseq.NT.samples))/2))

dataFilt <- dataFilt[filter,c(htseq.TP.samples, htseq.NT.samples)]

metadata <- colData(htseq.summarized.exp)
metadata <- metadata[colnames(dataFilt),]
metadata <- as.data.frame(metadata)
metadata.treament <- dplyr::select(metadata, c('treatments'))
metadata <- dplyr::select(metadata, -c('treatments'))
metadata <- metadata[metadata$sample_type != 'Recur',]

extension <- get_IDs(dataFilt)
extension <- dplyr::select(extension, -c(patient, sample))
metadata <- merge(metadata, extension, by='barcode')
rownames(metadata) <- metadata$barcode
metadata <- metadata[colnames(dataFilt),]
metadata$mirna_sample <- c(mirna.TP.samples, mirna.NT.samples)

# Error in plot.window(...) : need finite 'ylim' values; Error in plotting

# batch.corrected <- TCGAbatch_Correction(tabDF = dataFilt, batch.factor = 'Plate')

# Added plotting options; plot=FALSE to skip plotting

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

write.table(metadata, '../Data/TCGA-BRCA-COV.tsv', sep='\t', row.names=FALSE)

saveRDS(dataFilt, file='../Data/TCGA-BRCA-GENE-COUNTS.RDS')
