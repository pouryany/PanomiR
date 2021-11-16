#' Covariates/metadata of The Cancer Genome Atlas Liver Hepatocellular
#' Carcinoma gene expression data.
#'
#' This dataset is filtered to contain only sample with miRNA and mRNA gene
#' expression profiles.
#'
#' @format A data frame with 417 rows and 73 variables.
#' @source \url{https://doi.org/10.1016/j.cell.2017.05.046}
"COV_TCGA_LIHC"


#' The gene expression data The Cancer Genome Atlas Liver Hepatocellular
#' Carcinoma
#'
#' This dataset is filtered to contain only sample with miRNA and mRNA gene
#' expression profiles.
#'
#' @format A matrix with 14120 rows and 417 columns.
#' @source \url{https://doi.org/10.1016/j.cell.2017.05.046}
"TCGA_LIHC"




#' Canonical pathways from Molecular Signatures Database, MsigDb V6.2
#'
#' @format A list of 1329 pathways
#' @source \url{http://www.gsea-msigdb.org/gsea/index.jsp}
"msigdb_c2"


#' A table of gene-pathway association. Created based on the pathways of MSigDB.
#'
#' @format A matrix with 3 columns and 76926 rows:
#' \describe{
#'   \item{Pathway}{An MSigDB annotated pathway}
#'   \item{ENTREZID}{The ENTREZID of a gene belonging to the pathway}
#'   \item{ENSEMBL}{The ENSEMBL of a gene belonging to the pathway}
#' }
"path_gene_table"


#' The pathway correlation network of canonical pathways from the MSigDB.
#'
#' Each row represents two pathways, their overlap, their correlation
#' coefficient, and respective p-value of correlation.
#'
#' Please cite the appropriate publication if you use pcxn.
#'
#' https://doi.org/10.1371/journal.pcbi.1006042
#'
#' @format A matrix with 6 columns and 882456 rows:
#' \describe{
#'   \item{Pathway.A}{An MSigDB annotated pathway}
#'   \item{Pathway.B}{The ENTREZID of a gene belonging to the pathway}
#'   \item{Overlap.coefficient}{Jaccard index of the two pathways}
#'   \item{PathCor}{Correlation coefficient of the two pathways}
#'   \item{p.value}{correlation p.value}
#'   \item{p.Adjust}{FDR corrected p.value of correlation}
#' }
#' @source \url{https://pcxn.org}
"pcxn"


#' A processed list of miRNA target gene sets from the TargetScan dataset.
#' Each list item is a list of genes targeted by the respective miRNA family
#'
#' The interactions are filtered to only human interactions.
#'
#' The interactions are filtered to have a
#' Cumulative weighted context++ score of < 0
#' @format A list of 439 items
#' @source \url{http://www.targetscan.org/vert_72/}
"targetScan_00"


#' A processed list of miRNA target gene sets from the TargetScan dataset.
#' Each list item is a list of genes targeted by the respective miRNA family
#'
#' The interactions are filtered to only human interactions.
#'
#' The interactions are filtered to have a
#' Cumulative weighted context++ score of < -0.1
#' @format A list of 439 items
#' @source \url{http://www.targetscan.org/vert_72/}
"targetScan_01"


#' A processed list of miRNA target gene sets from the TargetScan dataset.
#' Each list item is a list of genes targeted by the respective miRNA family
#'
#' The interactions are filtered to only human interactions.
#'
#' The interactions are filtered to have a
#' Cumulative weighted context++ score of < -0.2
#' @format A list of 439 items
#' @source \url{http://www.targetscan.org/vert_72/}
"targetScan_02"

#' A processed list of miRNA target gene sets from the TargetScan dataset.
#' Each list item is a list of genes targeted by the respective miRNA family
#'
#' The interactions are filtered to only human interactions.
#'
#' The interactions are filtered to have a
#' Cumulative weighted context++ score of < -0.3
#' @format A list of 439 items
#' @source \url{http://www.targetscan.org/vert_72/}
"targetScan_03"
