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
#'   enrichment p-values
#' @export
miRNAPathwayEnrichment <- function(mirSets,
                                   pathwaySets,
                                   geneSelection = NULL,
                                   mirSelection = NULL,
                                   fromID = "ENSEMBL",
                                   toID = "ENTREZID",
                                   minPathSize = 9,
                                   numCores = 1,
                                   outDir = "",
                                   saveOutName = NULL) {
    if (substring(outDir, nchar(outDir)) != "/") {
        outDir <- paste0(outDir, "/")
    }
    if (!dir.exists(outDir)) {
        stop("Output directory does not exist.")
    }
    # select pathways with minimum set size
    pathsSel <- vapply(pathwaySets, length, numeric(1))
    pathwaySets <- pathwaySets[pathsSel > minPathSize]
    pathsRef <- Reduce(union, pathwaySets)

    # select miRNAs with targets in pathways
    mirSets <- lapply(mirSets, function(X) {
        X[X %in% pathsRef]
    })

    # select pathways with selected genes of interest
    if (!is.null(geneSelection)) {
        geneDF <- clusterProfiler::bitr(
            geneSelection,
            fromType = fromID,
            toType = toID,
            OrgDb = org.Hs.eg.db::org.Hs.eg.db
        )
        pathwaySets <- lapply(pathwaySets, function(X) {
            X[X %in% geneDF[, c(toID)]]
        })

        mirSets <- lapply(mirSets, function(X) {
            X[X %in% geneDF[, c(toID)]]
        })
        pathsRef <- Reduce(union, pathwaySets)
    }

    # select miRNAs of interest
    if (!is.null(mirSelection)) {
        mirSets <- mirSets[names(mirSets) %in% mirSelection]
    }

    selVec <- vapply(mirSets, length,numeric(1))
    mirSets <- mirSets[selVec > minPathSize]
    iterator <- (merge(names(mirSets), names(pathwaySets)))
    iterator <- iterator %>% dplyr::mutate_all(., as.character)
    all <- length(pathsRef)

    # find enrichment p-value of each miRNA target set and each pathway set
    enrichs <- parallel::mclapply(
        seq_len(nrow(iterator)),
        function(Y) {
            X <- iterator[Y, ]
            q <- length(intersect(
                unlist(pathwaySets[X[[2]]]),
                unlist(mirSets[X[[1]]])
            ))
            m <- length(unlist(mirSets[X[[1]]]))
            n <- all - m
            k <- length(unlist(pathwaySets[X[[2]]]))
            pval <- stats::phyper(
                q - 1, m, n, k, lower.tail = FALSE, log.p = FALSE
            )
            return(c(
                "pval" = pval,
                "Intersect" = q,
                "mirset_Size" = m,
                "not_mirset" = n,
                "pathway_Size" = k,
                "ratio_in" = q / (k - q),
                "ratio_out" = (m - q) / n,
                "ratio_ratios" = (q / (k - q)) / ((m - q) / n)
            ))
        },
        mc.cores = numCores
    )

    temp <- do.call(rbind, enrichs)
    temp <- as.data.frame(temp)
    iterator <- cbind(iterator, temp)

    if (!is.null(saveOutName)) {
        saveRDS(iterator, paste0(outDir, saveOutName))
    }
    return(iterator)
}
