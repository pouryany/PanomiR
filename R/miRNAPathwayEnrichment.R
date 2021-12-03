#' Enrichment Probability Of miRNAs
#'
#' Outputs enrichment probability of miRNAs based on pathway clusters.
#'
#' @param mirSets Table of miRNAs and a list of their interactions with
#'   genes in ENTREZ ID.
#' @param pathwaySets Table of pathways and a list of their interactions
#'   with genes in ENTREZ ID.
#' @param geneSelection Table of genes with dtype; if not NULL, select only
#'   genes from a given table.
#' @param mirSelection Table of miRNA names; if not NULL, select only miRNAs
#'   from given table.
#' @param fromID ID of genes in geneSelection.
#' @param toID ID of genes used in pcxn and pathways set.
#' @param minPathSize Filter out pathways with sets less than given value.
#' @param numCores Number of CPU cores to use, must be at least one.
#' @param outDir Output directory.
#' @param saveOutName If not NULL, saves output as RDS using save name.
#' @return Table of enrichment, each row contains mirna-pathway and its
#'   enrichment p-values.
#' @examples
#' data(msigdb_c2)
#' data(targetScan_03)
#' miRNAPathwayEnrichment(targetScan_03[1:20],msigdb_c2[1:20])
#' @export
miRNAPathwayEnrichment <- function(mirSets,
                                    pathwaySets,
                                    geneSelection = NULL,
                                    mirSelection = NULL,
                                    fromID = "ENSEMBL",
                                    toID = "ENTREZID",
                                    minPathSize = 9,
                                    numCores = 1,
                                    outDir = ".",
                                    saveOutName = NULL) {
    if (!dir.exists(outDir)) {
        stop("Output directory does not exist.")
    }

    # select pathways with minimum set size
    pathsSel <- vapply(pathwaySets, length, numeric(1))
    pathwaySets <- pathwaySets[pathsSel > minPathSize]
    pathsRef <- Reduce(union, pathwaySets)

    # select miRNAs with targets in pathways
    mirSets <- alignToUniverse(mirSets, pathsRef)
    # select pathways with selected genes of interest
    if (!is.null(geneSelection)) {
        geneDF <- clusterProfiler::bitr(geneSelection, fromType = fromID,
            toType = toID, OrgDb = org.Hs.eg.db::org.Hs.eg.db)

        pathwaySets <- alignToUniverse(pathwaySets, geneDF[, c(toID)])
        mirSets     <- alignToUniverse(mirSets, geneDF[, c(toID)])
        pathsRef    <- Reduce(union, pathwaySets)
    }
    # select miRNAs of interest
    if (!is.null(mirSelection)) {
        mirSets <- mirSets[names(mirSets) %in% mirSelection]
    }
    selVec   <- vapply(mirSets, length, numeric(1))
    mirSets  <- mirSets[selVec > minPathSize]

    iterator <- enrichAllPairs(mirSets, pathwaySets, pathsRef, numCores)

    if (!is.null(saveOutName)) {
        saveRDS(iterator, paste0(outDir, "/", saveOutName))
    }
    return(iterator)
}
