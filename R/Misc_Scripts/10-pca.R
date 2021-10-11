library(ggplot2)
library(ggfortify)

pathwaySummaryStats <- readRDS("../Data/TCGA-LIHC-x2-all-PathwaySummaryStats.RDS")
pathwaySummaryStats <- readRDS("../Data/TCGA-LIHC-x2-all-PathwaySummaryStats-NEW.RDS")
covariates <- read.csv('../Data/TCGA-LIHC-COV.csv', row.names = 1)

pathwaySummaryStats <- as.data.frame(t(pathwaySummaryStats))
pca_res <- prcomp(pathwaySummaryStats, scale. = TRUE)

autoplot(pca_res, data = covariates, colour = 'shortLetterCode')
