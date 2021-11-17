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
