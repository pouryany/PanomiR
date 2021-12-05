#' Canonical pathways from Molecular Signatures Database, MsigDb V6.2
#'
#' @format A list of 1329 pathways
#' @source \url{http://www.gsea-msigdb.org/gsea/index.jsp}
#' @usage data(msigdb_c2)
#' @examples
#' data(msigdb_c2)
"msigdb_c2"

#' A table of gene-pathway association. based on the pathways of MSigDB.
#'
#' @format A matrix with 3 columns and 76926 rows:
#' \describe{
#'   \item{Pathway}{An MSigDB annotated pathway}
#'   \item{ENTREZID}{The ENTREZID of a gene belonging to the pathway}
#'   \item{ENSEMBL}{The ENSEMBL of a gene belonging to the pathway}
#' }
#' @usage data(path_gene_table)
#' @examples
#' data(path_gene_table)
"path_gene_table"

#' A processed list of miRNA target gene sets from the TargetScan dataset.
#' Each list item is a list of genes targeted by the respective miRNA family
#'
#' The interactions are filtered to only human interactions.
#'
#' The interactions are filtered to have a
#' Cumulative weighted context++ score of < -0.3
#' @format A list of 439 items
#' @source \url{http://www.targetscan.org/vert_72/}
#' @usage data(targetScan_03)
#' @examples
#' data(targetScan_03)
"targetScan_03"

#' Readouts and datasets for minimal reproducible examples of the PanomiR.
#'
#' The item miniEnrich is a reduced representation of the TargetScan
#' For full table use miRNAPathwayEnrichment function in the package
#' along with msigdb_c2 and targetScan_03 datasets.
#'
#' These datasets include reduced representation of TCGA LIHC data
#' for reproducing the pipeline. doi: 10.1016/j.cell.2017.05.046
#'
#' A reduced representation of PCxN is provided. For full dataset and method
#' please refer to pcxn.org or https://doi.org/10.1371/journal.pcbi.1006042
#'
#' @format A list of 5:
#' \describe{
#'   \item{mini_LIHC_Exp}{a reduced expression dataset from TCGA LIHC data}
#'   \item{mini_LIHC_Cov}{a reduced covariates dataset from TCGA LIHC data}
#'   \item{miniEnrich}{a reduced table of miRNA-pathway enrichment, TargetScan.}
#'   \item{miniDEP}{Differentially activated pathways from reduced TCGA LIHC}
#'   \item{miniPCXN}{reduced representation of PCXN network}
#'   \item{miniPathClusts}{miniDEP mapped to miniPCXN}
#' }
#' @usage data(miniTestsPanomiR)
#' @examples
#' data(miniTestsPanomiR)
"miniTestsPanomiR"
