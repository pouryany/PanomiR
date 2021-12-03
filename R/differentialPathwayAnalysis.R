#' Differential Expression Analysis For Pathways
#'
#' Performs differential expression analysis for pathways using LIMMA package
#' with gene counts
#'
#' @param geneCounts Gene counts, rows refer to genes and columns to samples.
#' @param pathways Pathways table, containing pathway names and genes with id
#'   specified.
#' @param covariates Covariates/metadata file; rows matches the columns of
#'   geneCounts.
#' @param condition Condition to be examined (tumor vs normal etc); must exist
#'   in covariates column.
#' @param adjustCovars Adjustment covariates like batch; if NULL,
#'   no adjustments performed.
#' @param covariateCorrection If TRUE, performs covariates detection and
#'   correction; requires **adjustCovars**; (limma).
#' @param quantileNorm If TRUE, performs quantile normalization on pathway
#'   summary statistics; from *preprocess* package.
#' @param outDir Output directory.
#' @param saveOutName If not NULL, saves output as RDS using save name,
#'   if NULL, does not save output.
#' @param id ID matching genes to pathways; rownames of geneCounts.
#' @param deGenes If not NULL, add t-scores to pathways summary statistics;
#'   filter by genes t-scores.
#' @param minPathSize Minimum pathway size.
#' @param method Define method to use for pathway summary statistics;
#'   specifications in documentations.
#' @param trim Filter pathways with mean less than trim threshold
#'   in pathway summary statistics.
#' @param geneCountsLog If TRUE, log(geneCounts).
#' @param contrastConds Provide a contrast expression to be used in Limma
#'   comparison. This is necessary if you have more than two levels in the
#'   condition covariate.
#' @return List containing differentially expressed pathways as DEP and pathway
#'   summary statistics as pathwaySummaryStats.
#' @export
differentialPathwayAnalysis <- function(geneCounts,
                                        pathways,
                                        covariates,
                                        condition,
                                        adjustCovars = NULL,
                                        covariateCorrection = FALSE,
                                        quantileNorm = FALSE,
                                        outDir = ".",
                                        saveOutName = NULL,
                                        id = "ENSEMBL",
                                        deGenes = NULL,
                                        minPathSize = 10,
                                        method = "x2",
                                        trim = 0.025,
                                        geneCountsLog = TRUE,
                                        contrastConds = NA) {
    if (substring(outDir, nchar(outDir)) != "/") {
        outDir <- paste0(outDir, "/")
    }
    if (!dir.exists(outDir)) {
        stop("Output directory does not exist.")
    }

    # select pathways with genes in the gene count data and a minimum
    # pathway set size
    pathways <- as.data.frame(pathways)
    genesPathways <- pathways[pathways[, id] %in% rownames(geneCounts), ]

    genesPathways %<>% dplyr::group_by(., Pathway) %>%
        dplyr::summarise(., n = dplyr::n()) %>%
        dplyr::filter(., n >= minPathSize)

    pathways <- pathways %>%
        dplyr::filter(., Pathway %in% genesPathways$Pathway)

    # use de genes as a filter
    if (!is.null(deGenes)) {
        tScores <- deGenes %>%
            dplyr::mutate(., !!id := rownames(deGenes)) %>%
            dplyr::select(., c(ENSEMBL, t))
    }

    # log gene counts if gene counts are not log-transformed yet; essential for
    # path summary statistics
    if (geneCountsLog == TRUE) {
        geneCounts <- log(geneCounts)
    }

    # generate pathway summary statistics
    pathwaySummaryStats <-
        pathwaySummary(geneCounts, pathways, id = id, method = method,
            deGenes = deGenes, zNormalize = FALSE, trim = trim,
            tScores = tScores)

    # filter pathways with na values in the pathway summary statistics and
    # z-normalize the pathway summary statistics
    pathwaySummaryStats <-
        pathwaySummaryStats[rowSums(is.na(pathwaySummaryStats)) == 0, ]

    pathwaySummaryStats <- apply(pathwaySummaryStats, 2, function(x) {
        (x - mean(x)) / stats::sd(x)
    })

    # perform quantile normalization if needed
    # Add importing
    if (quantileNorm == TRUE) {
        pathwayNames <- rownames(pathwaySummaryStats)

        pathwaySummaryStats <-
            preprocessCore::normalize.quantiles(pathwaySummaryStats)

        rownames(pathwaySummaryStats) <- pathwayNames
        colnames(pathwaySummaryStats) <- rownames(covariates)
    }

    # set factors in covariates if needed

    # perform covariates correction, create design matrix for limma DE analysis;
    # if no adjustCovars are available, then design matrix only consider the
    # condition in question.
    fitResiduals <- NULL
    if (covariateCorrection == TRUE) {
        stop("Under development. Please use covariateCorrection = FALSE option")
    } else {
        if (is.null(adjustCovars)) {
            conditionsDat <- as.data.frame(covariates[, condition])
            colnames(conditionsDat) <- condition
            rownames(conditionsDat) <- rownames(covariates)
            designMat <- getDesignMatrix(conditionsDat, intercept = FALSE)
        } else {
            designMat <-
                getDesignMatrix(covariates[, c(condition, adjustCovars),
                                        drop = FALSE], intercept = FALSE)

            designMat$design <-
                designMat$design[, linColumnFinder(designMat$design)$indepCols]

            fitResiduals <- getResidual(covariates, adjustCovars,
                                    pathwaySummaryStats)
        }
    }

    conditionsTypes <- as.character(unique(covariates[, condition]))
    if (is.na(contrastConds)) {
        if (length(conditionsTypes) > 2) {
            stop("Please compare only 2 conditions at once.")
        }
        cond1 <- paste0(condition, conditionsTypes[1])
        cond2 <- paste0(condition, conditionsTypes[2])
        contrastsName <- paste0(cond1, "-", cond2)
    } else {
        contrastsName <- contrastConds
    }
    # limma DE analysis
    tT <- getDiffExpTable(pathwaySummaryStats, designMat, contrastsName)

    output <- list("DEP" = tT, "pathwaySummaryStats" = pathwaySummaryStats,
        "contrast" = contrastsName, "PathwayResiduals" = fitResiduals)

    if (!is.null(saveOutName)) {
        saveRDS(output, paste0(outDir, saveOutName))
    }
    return(output)
}
