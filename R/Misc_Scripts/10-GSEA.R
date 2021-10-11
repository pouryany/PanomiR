library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

# data(geneList)
# ego3 <- gseGO(geneList     = geneList,
#               OrgDb        = org.Hs.eg.db,
#               ont          = "CC",
#               nPerm        = 1000,
#               minGSSize    = 100,
#               maxGSSize    = 500,
#               pvalueCutoff = 0.05,
#               verbose      = FALSE)
# 
# gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
# c5 <- read.gmt(gmtfile)
# 
# egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)

geneList <- readRDS('../Data/TCGA-LIHC-DEG.RDS')
# geneList <- geneList[geneList$adj.P.Val < 0.1,]
geneList <- geneList[geneList$adj.P.Val < 0.05,]
genes <- rownames(geneList)
gene.ids <- bitr(genes, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
geneList <- geneList[gene.ids$ENSEMBL,]
geneList <- geneList[order(geneList$t, decreasing = T),]
geneList <- geneList$t
names(geneList) <- gene.ids$ENTREZID

pathways <- readRDS('../../Data/preprocessed/MSigDBPathGeneTab.RDS')
pathways <- pathways[,c(1,2)]

temp <- GSEA(geneList, 
             TERM2GENE = pathways, 
             verbose = F,
             pAdjustMethod = 'BH',
             pvalueCutoff = 1)

temp1 <- enricher(as.character(gene.ids$ENTREZID), TERM2GENE = pathways)
