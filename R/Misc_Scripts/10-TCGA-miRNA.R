rm(llist = ls())
library(TCGAbiolinks)
library(limma)
library(dplyr)
library(SummarizedExperiment)
library(miRBaseConverter)

# GDCdownload(query3, directory = '../GDCdata')
# mirna.data <- GDCprepare(query2, directory = '../GDCdata')
# mirna.names <- miRNA_NameToAccession(mirna.data$miRNA_ID)
# mirna.ids <- miRNA_AccessionToName(mirna.names$Accession, targetVersion = 'v21')


# Download from http://firebrowse.org/; miRseq_Mature_Preprocess (MD5)
mirna.data <- read.csv('../Data/LIHC_miRseq_mature_RPM.txt', sep='\t')
rownames(mirna.data)<- sapply(as.character(mirna.data$Gene), function(X){
  return(substr(X, 1, nchar(X)-13))
})
mirna.data <- mirna.data %>% dplyr::select(., -c('Gene'))
colnames(mirna.data) <- gsub('\\.', '-', colnames(mirna.data))
mirna.data[is.na(mirna.data)] <- 0

# 
# filter <- (rowSums(mirna.data[,mirna.TP.samples] > 1) >= floor(length(mirna.TP.samples)/2) &
#              rowSums(mirna.data[,mirna.NT.samples] > 1) >= floor(length(mirna.NT.samples)/2))

covariates <- read.csv('../Data/TCGA-LIHC-COV.csv', row.names = 1)
mirna.TP.samples <- as.character(covariates[covariates$shortLetterCode == 'TP',]$mirna_sample)
mirna.NT.samples <- as.character(covariates[covariates$shortLetterCode == 'NT',]$mirna_sample)

filter <- rowSums(mirna.data >= 1) >= floor(length(c(mirna.TP.samples, mirna.NT.samples))/20)

temp.mirna.data <- mirna.data[filter, c(mirna.TP.samples, mirna.NT.samples)]
# mirna.data <- mirna.data[,c(mirna.TP.samples, mirna.NT.samples)]

# saveRDS(mirna.data, '../Data/TCGA-LIHC-miRNA-COUNTS.RDS')

IDs <- get_IDs(temp.mirna.data)
type <- factor(c(rep('TP', length(mirna.TP.samples)),rep('NT',length(mirna.NT.samples))))
batch <- factor(IDs$plate)
design <- model.matrix(~0+type+batch)
fit <- lmFit(log(temp.mirna.data+3), design)
contr <- makeContrasts(typeTP-typeNT, levels = colnames(coef(fit)))
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit, trend=T)
tt <- topTable(fit, coef = 1, number=Inf)
saveRDS(tt, file=paste0('../Data/TCGA-LIHC-DEmiRNAs.RDS'))

# saveRDS(tt, file=paste0('../Data/TCGA-BRcA-DEmiRNAs.RDS'))


