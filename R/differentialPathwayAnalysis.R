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
#' @examples
#'
#' data("path_gene_table")
#' data("miniTestsPanomiR")
#'
#' differentialPathwayAnalysis(geneCounts = miniTestsPanomiR$mini_LIHC_Exp,
#' pathways =  path_gene_table,
#' covariates = miniTestsPanomiR$mini_LIHC_Cov,
#' condition = 'shortLetterCode')
#'
#' @export
differentialPathwayAnalysis <- function(geneCounts, pathways, covariates,
            condition, adjustCovars = NULL, covariateCorrection = FALSE,
            quantileNorm = FALSE, outDir = ".", saveOutName = NULL,
            id = "ENSEMBL", deGenes = NULL, minPathSize = 10, method = "x2",
            trim = 0.025, geneCountsLog = TRUE, contrastConds = NA) {

    outDir   <- .dirChecks(outDir)
    pathways <- .pathwayCleaner(pathways, id, geneCounts, minPathSize)
    tScores  <- .tScoreMaker(deGenes, id)

    if (geneCountsLog == TRUE) {
        geneCounts <- log(geneCounts)
    }

    pathSumStats <- pathwaySummary(geneCounts, pathways, id, zNormalize = FALSE,
        method = method,  deGenes = deGenes, trim = trim, tScores = tScores)
    pathSumStats <- .pathSumCleaner(pathSumStats)
    pathSumStats <- .qNormalizer(quantileNorm, pathSumStats, covariates)

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
            designMat <- getDesignMatrix(covariates[, c(condition,
                            adjustCovars), drop = FALSE], intercept = FALSE)
            designMat$design <-
                designMat$design[, linColumnFinder(designMat$design)$indepCols]
            fitResiduals <- getResidual(covariates, adjustCovars, pathSumStats)
        }
    }
    contrastsName <- .constrastHelper(covariates, condition, contrastConds)
    tT            <- getDiffExpTable(pathSumStats, designMat, contrastsName)
    output        <- list("DEP" = tT, "pathwaySummaryStats" = pathSumStats,
                "contrast" = contrastsName, "PathwayResiduals" = fitResiduals)

    if (!is.null(saveOutName)) {
        saveRDS(output, paste0(outDir, saveOutName))
    }
    return(output)
}
