#' Prioritize miRNA
#'
#' Outputs a table of miRNA (ordered) with respective p-values derived from
#'   method for prioritization
#'
#' @param enriches0 miRNA-pathway enrichment dataset obtained from
#'   miRNAPathwayEnrichment.
#' @param pathwayClusters Pathway clusters, obtained from
#'   MappingPathwaysClusters.
#' @param method Vector of methods (pCut, AggInv, AggLog, sumz, sumlog).
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
#' @return Table of miRNA and p-values, each row contains a miRNA and its
#'   associated p-values from the methods.
#' @export
prioritizeMicroRNA <- function(enriches0,
                               pathwayClusters,
                               method = "AggInv",
                               methodThresh = NULL,
                               enrichmentFDR = 0.25,
                               topClust = 2,
                               sampRate = 1000,
                               outDir = "",
                               dataDir = "",
                               saveSampling = TRUE,
                               runJackKnife = TRUE,
                               saveJackKnife = FALSE,
                               numCores = 1,
                               saveCSV = TRUE,
                               prefix = "") {
    if (substring(outDir, nchar(outDir)) != "/") {
        outDir <- paste0(outDir, "/")
    }
    if (!dir.exists(outDir)) {
        warning("Output directory does not exist.")
        dir.create(outDir, recursive = TRUE)
    }

    if (outDir == "/") {
        outDir <- ""
    }

    if (substring(dataDir, nchar(dataDir)) != "/") {
        dataDir <- paste0(dataDir, "/")
    }
    if (!dir.exists(dataDir)) {
        stop("Data directory does not exist.")
    }

    if (dataDir == "/") {
        dataDir <- ""
    }

    output <- list()

    # count miRNA-pathway enrichment with p-value less than threshold
    enriches <- enriches0 %>% dplyr::filter(., Intersect != 0)
    enriches %<>% dplyr::group_by(., y) %>%
        dplyr::mutate(., path_fdr = stats::p.adjust(pval, method = "fdr"))

    enriches <- enriches %>%
        dplyr::mutate(., hit = ifelse(path_fdr < enrichmentFDR, 1, 0))

    # Need to fix this function to be able to work on specific clusters
    # instead of all top clusters.
    for (clustNo in seq_len(topClust)) {
        clustName <- paste0("Cluster", clustNo)

        print(paste0("Working on ", clustName, "."))

        # select pathways in cluster
        pathways <- as.character(
            pathwayClusters[pathwayClusters$cluster == clustNo, ]$Pathway
        )
        nPaths <- length(pathways)

        # formulate number of miRNA-pathway enrichment with p-value less
        # than threshold for each miRNA
        tempEnrich <- enriches[enriches$y %in% pathways, ]
        selector <- tempEnrich %>%
            dplyr::group_by(x) %>%
            dplyr::summarise(., "cluster_hits" = sum(hit))

        # perform p-value aggregation based on methodlogy provided
        for (i in 1:length(method)) {
            m <- method[i]

            print(paste0("Performing ", m, " function."))

            fn <- get(paste0(m, "Fn"))
            coverFn <- get(paste0(m, "CoverFn"))

            if (!is.null(methodThresh)) {
                mThresh <- methodThresh[i]
                temp <- fn(
                    enriches = enriches0,
                    pathways,
                    isSelector = TRUE,
                    thresh = mThresh
                )
            } else {
                temp <- fn(enriches = enriches0, pathways, isSelector = TRUE)
            }

            mSelector <- temp$selector
            mEnriches0 <- temp$enriches0

            if (nrow(mSelector) < 3) {
                print(paste0(
                    "Skipping ",
                    m,
                    " function due to low number of miRNA after filter"
                ))
                next
            }

            enrichNull <- mEnriches0 %>% dplyr::filter(., x %in% mSelector$x)

            samplingDataDir <- paste0(dataDir, prefix, "Sampling_Data/")

            if (saveSampling == TRUE) {
                if (!dir.exists(samplingDataDir)) {
                    dir.create(samplingDataDir, recursive = TRUE)
                }
            }

            samplingDataFilename <- paste0(
                prefix,
                m,
                "_",
                sampRate,
                "_samples.RDS"
            )

            samplingDataFile <- paste0(
                samplingDataDir,
                "/",
                samplingDataFilename
            )

            if (m %in% c("aggInv", "aggLog")) {
                # perform sampling
                samplingData <- samplingDataBase(enrichNull,
                                                 mSelector,
                                                 sampRate,
                                                 fn,
                                                 nPaths,
                                                 samplingDataFile,
                                                 jackKnife = FALSE,
                                                 saveSampling = saveSampling,
                                                 numCores = 8
                )

                mSelector <- methodProbBase(
                    samplingData = samplingData[[paste0("SampSize_", 100)]],
                    selector = mSelector,
                    m = m,
                    nPaths = nPaths,
                    coverFn = coverFn
                )
            } else {
                names(mSelector)[2:3] <- paste0(m, "_", names(mSelector)[2:3])
            }

            print(paste0(m, " Method Done"))

            # perform jack-knife
            if ((runJackKnife == TRUE) && (m != "sumz") && (m != "sumlog")) {
                samplingData <- samplingDataBase(enrichNull,
                                                 mSelector,
                                                 sampRate,
                                                 fn,
                                                 nPaths,
                                                 samplingDataFile,
                                                 jackKnife = FALSE,
                                                 saveSampling = saveSampling,
                                                 numCores = 8
                )
                mSelector <- jackKnifeBase(
                    selector = mSelector,
                    pathways = pathways,
                    enrichNull = enrichNull,
                    fn = fn,
                    jackKnifeData = samplingData[[paste0("SampSize_", 100)]],
                    m = m,
                    numCores = 8
                )

                print(paste0(m, " JackKnifing Method Done!"))

                mSelector <- mSelector[, c(1, 3, 2, 4)]
            } else {
                mSelector <- mSelector[, c(1, 3, 2)]
            }

            selector <- merge(selector, mSelector, all = TRUE)

            met <- paste0(m, "_pval")
            selector <- selector %>% dplyr::arrange(., !!rlang::sym(met))
            met2 <- paste0(m, "_fdr")
            selector <- selector %>% dplyr::mutate(
                .,
                !!met2 := stats::p.adjust(p = !!rlang::sym(met), method = "fdr")
            )
        }

        if (saveCSV == TRUE) {
            saveName <- paste0(
                prefix,
                sampRate,
                "_samples_clustNo_",
                clustNo, ".csv"
            )

            print(paste0(saveName, " saved!"))

            utils::write.csv(selector, paste0(outDir, saveName))
        }

        output[[clustName]] <- selector
    }
    return(output)
}
