#' Prioritize miRNA
#'
#' Outputs a table of miRNA ordered with respective p-values derived from
#'   method for prioritization
#'
#' @param enriches0 miRNA-pathway enrichment dataset obtained from
#'   miRNAPathwayEnrichment.
#' @param pathClust Pathway clusters, obtained from
#'   MappingPathwaysClusters.
#' @param method Vector of methods pCut, AggInv, AggLog, sumz, sumlog.
#' @param methodThresh Vector of methods threshold for each method in method,
#'   if NULL use default thresh values in method.
#' @param enrichmentFDR FDR cut-off calculating miRNA-pathway hits
#'   in the input cluster based on significant enrichment readouts.
#' @param topClust Top x clusters to perform miRNA prioritization on.
#' @param sampRate Sampling rate for CLT.
#' @param numCores Number of CPU cores to use, must be at least one.
#' @param outDir Output directory.
#' @param dataDir Data directory.
#' @param saveCSV If TRUE, saves CSV file for each cluster in topClust in
#'   outDir.
#' @param saveSampling If TRUE, saves sampling data as RDS for each cluster in
#'   topClust in dataDir.
#' @param runJackKnife If TRUE, jacknifing will be performed.
#' @param saveJackKnife If TRUE, saves jack-knifed sampling data as RDS for each
#'   cluster in topClust in dataDir.
#' @param prefix Prefix for all saved data.
#' @param autoSeed random permutations are generated based on predetermined
#'   seeds. TRUE will give identical results in different runs.
#' @return Table of miRNA and p-values, each row contains a miRNA and its
#'   associated p-values from the methods.
#' @examples
#' data("miniTestsPanomiR")
#'
#' prioritizeMicroRNA(enriches0 = miniTestsPanomiR$miniEnrich,
#'    pathClust = miniTestsPanomiR$miniPathClusts$Clustering,
#'    topClust = 1,
#'    sampRate = 50,
#'    method = c("aggInv"),
#'    saveSampling = FALSE,
#'    runJackKnife = FALSE,
#'    numCores = 1,
#'    saveCSV = FALSE)
#' @export
prioritizeMicroRNA <- function(enriches0, pathClust, method = "AggInv",
            methodThresh = NULL, enrichmentFDR = 0.25, topClust = 2,
            sampRate = 1000, outDir = ".", dataDir = ".", saveSampling = TRUE,
            runJackKnife = TRUE, saveJackKnife = FALSE, numCores = 1,
            saveCSV = TRUE, prefix = "", autoSeed = TRUE) {
    .checkAddressDirs(outDir, saveCSV, dataDir, saveSampling)
    output   <- list()
    enriches <- .cleanEnrichInput(enriches0, enrichmentFDR)

    for (clustNo in seq_len(topClust)) {
        clustName <- paste0("Cluster", clustNo)
        print(paste0("Working on ", clustName, "."))

        pathways <- as.character(pathClust[pathClust$cluster == clustNo,
                                    ]$Pathway)
        nPaths <- length(pathways)
        selector <- .makeSelector(enriches, pathways)

        for (i in seq_along(method)) {
            m <- method[i]
            print(paste0("Performing ", m, " function."))
            fn <- get(paste0(m, "Fn"))
            coverFn <- get(paste0(m, "CoverFn"))

            if (!is.null(methodThresh)) {
                mThresh <- methodThresh[i]
                temp <- fn(enriches = enriches0, pathways, isSelector = TRUE,
                    thresh = mThresh)
            } else {
                temp <- fn(enriches = enriches0, pathways, isSelector = TRUE)
            }
            mSelector   <- temp$selector
            mEnriches0  <- temp$enriches0
            if (nrow(mSelector) < 3) {
                print(paste0("Skipping ", m, " function:few miRNAs"))
                next
            }
            enrichNull <- mEnriches0 %>% dplyr::filter(., x %in% mSelector$x)
            sampDatFile <- .mkSampleDatDir(dataDir, prefix, saveSampling,
                                        m, sampRate)
            mSelector   <- .mSelectorMaker(m, enrichNull, mSelector, sampRate,
                        fn, nPaths, sampDatFile, saveSampling, numCores,
                        autoSeed, coverFn, runJackKnife, pathways)
            selector <- .expandSelector(selector, mSelector, m)
        }
        .saveCsvWriter(saveCSV, prefix, sampRate, clustNo, selector, outDir)
        output[[clustName]] <- selector
    }
    return(output)
}
