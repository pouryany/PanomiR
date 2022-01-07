utils::globalVariables(c(
    ".", ":=", "ENSEMBL", "ES", "ES2", "Intersect", "Pathway", "hit",
    "hit2", "k", "n", "path_fdr", "pval", "x", "y"
))

#' Pathway-Gene Associations
#'
#' Generates a table of pathways and genes associations.
#'
#' @param pathAdress Address to an RDS file containing list of pathways where
#'   each element is a list of genes similar to GMT format.
#' @param pathwayList If you wish to use a list of pathways instead of a file
#'   use this argument instead. The list must contain no NA values.
#' @param fromType gene annotation type used in your input data.
#' @param toType gene annotation type to be produced in the output.
#' @param outDir Address to save an RDS for a table of pathway-gene association
#' @return pathExpTab Table of pathway-gene association.
#' @examples
#' pathway1 <- c("125", "3099", "126")
#' pathway2 <- c("5232", "5230", "5162")
#' pathList <- list("Path1" = pathway1, "Path2" = pathway2)
#' res <- pathwayGeneTab(pathwayList = pathList)
#'
#' data(msigdb_c2)
#' pathwayGeneTab(pathwayList = msigdb_c2[1:2])
#' @export
pathwayGeneTab <- function(pathAdress = NA, pathwayList = NA,
                    fromType = "ENTREZID", toType = "ENSEMBL", outDir = NA) {
    if (!is.na(pathAdress)) {
        pathList <- readRDS(pathAdress)
    }

    if (is.na(pathAdress) && !any(is.na(pathwayList))) {
        pathList <- pathwayList
    }

    if (!is.na(pathAdress) && !any(is.na(pathwayList))) {
        stop("provide a valid input list.")
    }

    pathList2 <- lapply(
        pathList,
        function(x) {
            clusterProfiler::bitr(x, fromType = fromType, toType = toType,
                OrgDb = org.Hs.eg.db::org.Hs.eg.db)
        }
    )

    temp <- lapply(
        names(pathList2),
        function(x) {
            data.frame(
                Pathway = (x),
                pathList2[[x]]
            )
        }
    )
    temp <- do.call(rbind, (temp))
    pathExpTab <- tibble::as_tibble(temp)

    if (!is.na(outDir)) {
        saveRDS(pathExpTab, outDir)
    }

    return(pathExpTab)
}

#' Pathway-Gene Associations from GeneSet collections
#'
#' This function enables to utilize MSigDB packages and GSEABase objects
#'   to incorporate customized genesets into PanomiR.
#'
#' @param gsCollection An GSEABase gene set collection object
#' @param fromType gene annotation type used in your input data
#' @param toType gene annotation type to be produced in the output
#' @return A table of pathway-gene associations
#' @examples
#' data(gscExample)
#' tableFromGSC(gscExample)
#' @export
tableFromGSC <- function(gsCollection, fromType = "ENTREZID",
                                toType = "ENSEMBL") {
    gsList <- GSEABase::geneIds(gsCollection)
    pathwayGeneTab(pathwayList = gsList, fromType = fromType, toType = toType)
}

#' Pathway Summary Statistics
#'
#' Generates a table of pathway activity profiles per sample
#'
#' @param exprsMat Gene expression matrix with row names as genes and samples
#'   as columns.
#' @param pathwayRef Table of pathway-gene associations. Created from
#'   \code{\link{pathwayGeneTab}} function.
#' @param id Gene annotation type in the row name of gene expression data.
#' @param zNormalize Normalization of pathway summary score.
#' @param method Choice of how to summarize gene ranks into pathway statistics.
#' @param deGenes List of differentially expressed genes along with t-scores.
#'   Only necessary if working on Top 50\% summary method.
#' @param trim Percentage of top and bottom ranked genes to be excluded from
#'   pathway summary statistics.
#' @param tScores Argument for-top-50-percent-genes method.
#' @return pathExp Table of pathway activity profiles per sample.
#' @examples
#' pathTab <- tibble::tribble(
#'  ~Pathway, ~ENTREZID,  ~ENSEMBL,
#'  "Path1",   "125",      "ENSG00000196616",
#'  "Path1",   "3099",     "ENSG00000159399",
#'  "Path2",   "5230",     "ENSG00000102144",
#'  "Path2",   "5162",     "ENSG00000168291"
#'  )
#' exprsMat <- matrix(2 * (seq_len(12)), 4, 3)
#' rownames(exprsMat) <- pathTab$ENSEMBL
#' colnames(exprsMat) <- LETTERS[seq_len(3)]
#' pathwaySummary(exprsMat, pathTab, method = "x2")
#' @export

pathwaySummary <- function(exprsMat, pathwayRef, id = "ENSEMBL",
                    zNormalize = FALSE, method = FALSE, deGenes = NULL,
                    trim = 0, tScores = NULL) {
    # Current version only works with ENSEMBL.
    exprsMat <- tibble::rownames_to_column(as.data.frame(exprsMat), var = id)
    if (!is.null(deGenes)) {
        if (is.null(tScores)) {
            stop("Provide tscores/pvalues")
        }
        pathwayRef <- dplyr::inner_join(pathwayRef, tScores, by = c("ENSEMBL"))
        pathwayRef <- pathwayRef |> dplyr::group_by(Pathway) |>
            dplyr::filter(abs(t) >= stats::median(abs(t))) |>
            dplyr::select(-t)
    }
    if (method == "none") {
        exprsMat <- exprsMat |>
            dplyr::mutate_if(is.numeric, function(x) {
                x
                })
    } else if (method == "x") {
        exprsMat <- exprsMat |>
            dplyr::mutate_if(is.numeric, rank) |>
            dplyr::mutate_if(is.numeric, function(x) {
                x
                })
    } else if (method == "x2") {
        exprsMat <- exprsMat |>
            dplyr::mutate_if(is.numeric, rank) |>
            dplyr::mutate_if(is.numeric, function(x) {
                x * x
            })
    } else {
        stop("invalid choice of summarization function")
    }
    pathExpTab <- dplyr::inner_join(pathwayRef, exprsMat, by = id)
    pathExpTab <- pathExpTab |>
        dplyr::group_by(Pathway) |>
        dplyr::summarise_if(is.numeric, mean, na.rm = TRUE, trim = trim)

    pathExp <- as.data.frame(pathExpTab[, -1])
    rownames(pathExp) <- as.character(dplyr::pull(pathExpTab[, 1]))

    if (zNormalize) {
        pathExp <- apply(pathExp, 2, function(x) {
            (x - mean(x)) / stats::sd(x)
        })
    }
    return(pathExp)
}


### ==========================================================
#-----
#-----
#----- The following are miRNA-Pathway pval aggregations
#-----
#-----
### ==========================================================



#' Score miRNAs In a Cluster Of Pathways
#'
#' The function to count the number of enriched pathways for each miRNA.
#'
#' @param enriches Table of miRNA pathway enrichments.
#' @param pathways Queried pathways, e.g. cluster pathways.
#' @param isSelector Internal argument.
#' @param thresh Threshold from p-value cut-off.
#' @return P-value based scoring of miRNAs in a cluster of pathways.
#' @keywords internal
pCutFn <- function(enriches, pathways, isSelector, thresh = 0.05) {
    if (isSelector == TRUE) {
        enriches <- enriches |>
            dplyr::mutate(hit2 = ifelse(pval < thresh, 1, 0))
    }

    tempEnrich <- enriches[enriches$y %in% pathways, ]
    selector <- tempEnrich |>
        dplyr::group_by(x) |>
        dplyr::summarise(n = dplyr::n(), k = sum(hit2)) |>
        dplyr::arrange(x)

    if (isSelector == TRUE) {
        selector <- selector |> dplyr::filter(k > thresh * length(pathways))
        return(list("selector" = selector, "enriches0" = enriches))
    } else {
        return(selector)
    }
}


#' The function calculate targeting score of miRNA w.r.t to a cluster
#' of pathways via inverse normal method
#' @param enriches a table of miRNA pathway enrichments. Universe
#' @param pathways queried pathways. e.g. cluster pathways
#' @param isSelector internal argument
#' @param thresh internal argument
#' @return a  scoring of miRNAs in a cluster of pathways
#' @keywords internal
aggInvFn <- function(enriches, pathways, isSelector = TRUE, thresh = NULL) {
    if (isSelector == TRUE) {
        enriches <- enriches |> dplyr::mutate(ES2 = stats::qnorm(1 - pval))
        minES <- min(enriches$ES2[!is.infinite(enriches$ES2)])
        enriches <- enriches |>
            dplyr::mutate(ES2 = ifelse(is.infinite(enriches$ES2), minES,
                enriches$ES2))
    }

    tempEnrich <- enriches[enriches$y %in% pathways, ]
    selector <- tempEnrich |>
        dplyr::group_by(x) |>
        dplyr::summarise(n = dplyr::n(), k = mean(ES2)) |>
        dplyr::arrange(x)

    if (isSelector == TRUE) {
        return(list("selector" = selector, "enriches0" = enriches))
    } else {
        return(selector)
    }
}


#' The function calculate targeting score of miRNA w.r.t to a cluster
#' of pathways via log aggregation method.
#' @param enriches a table of miRNA pathway enrichments. Universe
#' @param pathways queried pathways. e.g. cluster pathways
#' @param isSelector internal argument
#' @param thresh internal argument
#' @return a  scoring of miRNAs in a cluster of pathways
#' @keywords internal
aggLogFn <- function(enriches, pathways, isSelector, thresh = 0) {
    enriches <- enriches |> dplyr::mutate(ES = -log(pval))
    if (isSelector == TRUE) {
        enriches <- enriches |> dplyr::mutate(ES = -log(pval))
    }

    tempEnrich <- enriches[enriches$y %in% pathways, ]
    selector <- tempEnrich |>
        dplyr::group_by(x) |>
        dplyr::summarise(n = dplyr::n(), k = mean(ES))

    if (isSelector == TRUE) {
        selector <- selector |>
            dplyr::filter(k * n >= thresh * length(pathways))
        return(list("selector" = selector, "enriches0" = enriches))
    } else {
        return(selector)
    }
}


#' The function calculate targeting score of miRNA w.r.t to a cluster
#' of pathways via sumz aggregation method.
#' @param enriches a table of miRNA pathway enrichments. Universe
#' @param pathways queried pathways. e.g. cluster pathways
#' @param isSelector internal argument
#' @param thresh internal argument
#' @return a  scoring of miRNAs in a cluster of pathways
#' @keywords internal
sumzFn <- function(enriches, pathways, isSelector, thresh = NULL) {
    enriches1 <- enriches |>
        dplyr::mutate(pval = ifelse(pval >= 0.999, 0.999, pval))

    enriches1 <- enriches1 |>
        dplyr::mutate(pval = ifelse(pval <= 1.0e-16, 1.0e-16, pval))

    tempEnrich <- enriches1[enriches1$y %in% pathways, ]
    aggPTab <- vector()
    for (i in unique(tempEnrich$x)) {
        temp <- tempEnrich[tempEnrich$x == i, ]
        tPVal <- metap::sumz(temp$pval)
        aggPTab <- rbind(aggPTab, c(i, tPVal$p, nrow(temp)))
    }

    selector <- tibble::tibble(
        "x" = aggPTab[, 1],
        "pval" = signif(as.numeric(aggPTab[, 2]), 4),
        "n" =  aggPTab[, 3]
    )

    if (isSelector == TRUE) {
        return(list("selector" = selector, "enriches0" = enriches))
    } else {
        return(selector)
    }
}


#' The function calculate targeting score of miRNA w.r.t to a cluster
#' of pathways via sumlog aggregation method.
#' @param enriches a table of miRNA pathway enrichments. Universe
#' @param pathways queried pathways. e.g. cluster pathways
#' @param isSelector internal argument
#' @param thresh internal argument
#' @return a  scoring of miRNAs in a cluster of pathways
#' @keywords internal
sumlogFn <- function(enriches, pathways, isSelector, thresh = NULL) {
    enriches1 <- enriches |> dplyr::mutate(pval = ifelse(pval >= 0.999,
                                    0.999, pval))
    enriches1 <- enriches1 |> dplyr::mutate(pval = ifelse(pval <= 1.0e-16,
                                    1.0e-16, pval))

    tempEnrich <- enriches1[enriches1$y %in% pathways, ]
    aggPTab <- vector()
    for (i in unique(tempEnrich$x)) {
        temp <- tempEnrich[tempEnrich$x == i, ]
        tPVal <- metap::sumlog(temp$pval)
        aggPTab <- rbind(aggPTab, c(i, tPVal$p, nrow(temp)))
    }

    selector <- tibble::tibble(
        "x" = aggPTab[, 1],
        "pval" = signif(as.numeric(aggPTab[, 2]), 4),
        "n" =  aggPTab[, 3]
    )

    if (isSelector == TRUE) {
        return(list("selector" = selector, "enriches0" = enriches))
    } else {
        return(selector)
    }
}


#' Internal function for modification of prioritization.
#' @param selector a prioritzation table
#' @param coverName a new column name
#' @return an updated scoring of miRNAs in a cluster of pathways
#' @keywords internal
pCutCoverFn <- function(selector, coverName) {
    selector <- selector |>
        dplyr::mutate(!!coverName := k / n)
    return(selector)
}


#' Internal function for modification of prioritization.
#' @param selector a prioritzation table
#' @param coverName a new column name
#' @return an updated scoring of miRNAs in a cluster of pathways
#' @keywords internal
aggInvCoverFn <- function(selector, coverName) {
    selector <- selector |>
        dplyr::mutate(!!coverName := k)
    return(selector)
}

#' Internal function for modification of prioritization.
#' @param selector a prioritzation table
#' @param coverName a new column name
#' @return an updated scoring of miRNAs in a cluster of pathways
#' @keywords internal
aggLogCoverFn <- aggInvCoverFn

#' Internal function for modification of prioritization.
#' @param selector a prioritzation table
#' @param coverName a new column name
#' @return an updated scoring of miRNAs in a cluster of pathways
#' @keywords internal
sumzCoverFn <- aggInvCoverFn

#' Internal function for modification of prioritization.
#' @param selector a prioritzation table
#' @param coverName a new column name
#' @return an updated scoring of miRNAs in a cluster of pathways
#' @keywords internal
sumlogCoverFn <- aggInvCoverFn


#### Working

#' Outputs a table of sampling data(rows are miRNA and cols are samples)
#'
#' @param enrichNull Enrichment dataset with x (miRNA), y (pathway) and pval
#'   (probability of observing x in pathway cluster).
#' @param selector Table with x(miRNA) in pathway cluster.
#' @param sampRate Sampling rate.
#' @param fn Methodology function.
#' @param nPaths Number of pathways in pathway cluster.
#' @param samplingDataFile If file exists, load. Else, perform random sampling
#' @param saveSampling If TRUE, data is saved.
#' @param jackKnife If TRUE, conduct sampling with one less pathway, used for
#'   jack knifing
#' @param numCores number of cores used
#' @param autoSeed random permutations are generated based on predetermined
#'   seeds. TRUE will give identical results in different runs.
#' @return Outputs of sampling data.
samplingDataBase <- function(enrichNull, selector, sampRate, fn, nPaths,
                        samplingDataFile, jackKnife = FALSE, saveSampling,
                        numCores = 1, autoSeed = TRUE) {
    if (!all(utils::hasName(selector, c("x")))) {
        stop("The selector table needs a column x (miRNA name)")
    }
    if (!all(utils::hasName(enrichNull, c("x", "y", "pval")))) {
        stop("Bad formatting in the enrichment table")
    }
    if (!file.exists(samplingDataFile)) {
        allPaths <- unique(enrichNull$y)
        if (jackKnife == TRUE) {
            nPathsTemp <- nPaths - 1
        }
        sampSizeVec <- c(nPaths, nPaths - 1, 100, 50)
        outList <- list()
        for (nPathsTemp in sampSizeVec) {
            temp <- parallel::mclapply(seq_len(sampRate), function(y) {
                if (autoSeed) {
                    withr::with_seed(y, {
                        nullPaths <- sample(allPaths, nPathsTemp,
                                        replace = FALSE)})
                }
                selNull <- fn(enriches = enrichNull, pathways = nullPaths,
                    isSelector = FALSE)
                return(selNull$k)
            }, mc.cores = numCores)
            # build null distribution of K
            temp <- do.call(rbind, temp)
            temp <- t(temp)
            rownames(temp) <- selector$x
            colnames(temp) <- vapply(seq_len(sampRate), function(y) {
                paste0("sample_", y)
            }, FUN.VALUE = "character")
            sampTag <- paste0("SampSize_", nPathsTemp)
            outList[[sampTag]] <- temp
        }
        if (saveSampling == TRUE) {
            saveRDS(outList, file = samplingDataFile)
            sayThis <- paste0(samplingDataFile, " saved.")
            message(sayThis)
        }
    } else {
        outList <- readRDS(samplingDataFile)
        sayThis <- paste0("Skipping sampling, ", samplingDataFile, " exists. ",
                        samplingDataFile, " loaded.")
        message(sayThis)
    }
    return(outList)
}

#' Outputs a table with col x, miRNA, probability of observing k
#'   against a random distribution of the cover of methodology
#'
#' @param samplingData Random distribution data.
#' @param selector Table with x(miRNA) in pathway cluster and observed
#'   k (depending on methodology).
#' @param m Method name.
#' @param nPaths Number of pathways used to generate the samplingData at
#'   each iteration. Default is set at 100.
#' @param coverFn Cover of methodology function.
#' @return Outputs a new selector table with col x, pval and cover.
methodProbBase <- function(samplingData,
                            selector,
                            m,
                            nPaths = 100,
                            coverFn = NULL) {
    if (!all(utils::hasName(selector, c("x", "k")))) {
        stop(
            "The selector needs a column x (miRNA) and a column k (miRNA hits)"
            )
    }

    # obtain means and sds for distribution, assume CLT
    means <- rowMeans(samplingData)
    sds <- apply(samplingData, 1, stats::sd)
    sds <- sds * 10 / sqrt(nPaths)

    pvalName <- paste0(m, "_pval")
    coverName <- paste0(m, "_cover")

    # obtain p-vals
    pVals <-
        stats::pnorm(selector$k, mean = means, sd = sds, lower.tail = FALSE)
    selector <- selector |>
        dplyr::mutate(!!pvalName := pVals) |>
        coverFn(coverName) |>
        dplyr::select(-c(k, n))
    return(selector)
}


#' Outputs a table with col x (miRNA), probability of observing
#' k (depending on methodology) against a random distribution with
#' jack-knifing of the pathway cluster (removing a pathway at a time)
#'
#' @param selector Table with x(miRNA) in pathway cluster and observed k
#'   (depending on methodology).
#' @param pathways Pathways in pathway cluster.
#' @param enrichNull Enrichment dataset with x (miRNA), y (pathway) and pval
#'   (probability of observing x in pathway cluster).
#' @param fn Methodology function.
#' @param jackKnifeData Random distribution data with jack-knifing
#'   (i.e. one less pathway)
#' @param m method name
#' @param numCores number of cores
#' @return Outputs a new selector table with col x, pval_jk
jackKnifeBase <- function(selector,
                            pathways,
                            enrichNull,
                            fn,
                            jackKnifeData,
                            m,
                            numCores = 1) {
    # obtain means and sds for distribution, assume CLT
    nPaths <- length(pathways)
    sampleMeans <- rowMeans(jackKnifeData)
    sampleSDs <- apply(jackKnifeData, 1, stats::sd)
    sampleSDs <- sampleSDs * 10 / sqrt(nPaths - 1)

    # remove one pathway at a time and obtain K for each miRNA
    temp1 <- parallel::mclapply(seq_along(pathways), function(x) {
        tempPathways <- pathways[-x]
        tempSelector <- fn(
            enriches = enrichNull, pathways = tempPathways,
            isSelector = FALSE
        )
        # obtain p-values using the means and sds obtain above
        pVals  <- stats::pnorm(tempSelector$k, mean = sampleMeans,
                        sd = sampleSDs, lower.tail = FALSE)

        return(pVals)
    }, mc.cores = numCores)

    # rows <- number of pathways; col <- number of miRNAs
    # contain p-values
    temp1 <- do.call(rbind, temp1)
    temp1 <- t(temp1)

    ## added from here
    jackKnifeName <- paste0(m, "_pval_jk")

    # obtain aggregate p-values
    means <- rowMeans(temp1)
    selector <- selector |>
        dplyr::mutate(!!jackKnifeName := means)

    return(selector)
}

#' Obtain Design Matrix
#'
#' Modified from covariates pipeline of Menachem Former. Imported from
#' \url{https://github.com/th1vairam/CovariateAnalysis}
#'
#' @param covariatesDataFrame Dataframe of covariates.
#' @param intercept intercept in the linear model.
#' @param reLevels TBA.
#' @return List containing a design matrix.
#' @examples
#' data(iris)
#' getDesignMatrix(iris)
#' @export
getDesignMatrix <- function(covariatesDataFrame, intercept = TRUE,
                            reLevels = list()) {
    covDF <- covariatesDataFrame
    rowNamesTemp <- rownames(covDF)
    colNamesTemp <- colnames(covDF)
    factorCovariateNames <- names(covDF)[vapply(covDF, is.factor, logical(1))]
    factorCovariateNames <- setdiff(factorCovariateNames,
        factorCovariateNames[!(factorCovariateNames %in% colnames(covDF))])

    numericCovariateNames <- setdiff(colNamesTemp, factorCovariateNames)

    # Ensure the factors are factor, and the quantitatives are numeric:
    covDF <- as.data.frame(lapply(colnames(covDF),
        function(column) {
            if (column %in% factorCovariateNames) {
                fac <- as.factor(covDF[, column])
                if (column %in% names(reLevels)) {
                    fac <- stats::relevel(fac, ref = reLevels[[column]])
                }
                return(fac)
            } else {
                return(as.numeric(covDF[, column]))
            }
        }))
    rownames(covDF) <- rowNamesTemp
    colnames(covDF) <- colNamesTemp

    tempHelper <- .getDesignMatHelper(factorCovariateNames, covDF, intercept)
    contra   <- tempHelper$contra
    covDF    <- tempHelper$covDF
    # Inspired by http://stackoverflow.com/questions/5616210/
    currentNAAction <- getOption("na.action")
    # Model matrix will now include "NA":
    options(na.action = "na.pass")

    if (intercept) {
        design <- stats::model.matrix(~., data = covDF, contrasts.arg = contra)
    } else {
        design <- stats::model.matrix(~ 0 + ., data = covDF,
                            contrasts.arg = contra)
    }

    rownames(design) <- rownames(covDF)
    options(na.action = currentNAAction)

    return(list(design = design, covariates = colNamesTemp,
        factorsLevels = lapply(contra, colnames),
        numericCovars = numericCovariateNames, covariatesDataFrame = covDF))
}

.getDesignMatHelper <- function(factorCovariateNames, covDF,
                                intercept) {
    contra <- NULL
    maxNumCat <- Inf
    catData <- covDF[, factorCovariateNames, drop = FALSE]
    if (ncol(catData) > 0) {
        numCats <- vapply(colnames(catData),
            function(col) nlevels(factor(catData[, col])), numeric(1))
        excludeCategoricalCols <-
            names(numCats)[numCats <= 1 | numCats > maxNumCat]

        if (!is.null(excludeCategoricalCols) &&
            length(excludeCategoricalCols) > 0) {
            warning(paste("Excluding categorical variables with less than 2",
                ifelse(is.infinite(maxNumCat), "",
                    paste(" or more than ", maxNumCat, sep = "")),
                " categories: ",
                paste(paste("'", excludeCategoricalCols,
                "'", sep = ""), collapse = ", "), sep = ""))

            factorCovariateNames <- setdiff(factorCovariateNames,
                excludeCategoricalCols)
            covDF <- covDF[, !(colnames(covDF) %in% excludeCategoricalCols),
                drop = FALSE]
        }
        # Inspired by http://stackoverflow.com/questions/4560459/
        # And, already ensured above that covDF[, factorCovariateNames] :
        # 1) fac is of type factor.
        # 2) fac is releveled as designated in reLevels.
        if (intercept) {
            contra <- lapply(
                factorCovariateNames,
                function(column) {
                    fac <- covDF[, column]
                    fac <- stats::contrasts(fac)
                }
            )
        } else {
            contra <- lapply(factorCovariateNames,
                function(column) {
                    fac <- covDF[, column]
                    fac <- stats::contrasts(fac, contrasts = FALSE)
                })
        }
        names(contra) <- factorCovariateNames
    }
    return <- list("contra" = contra, "covDF" = covDF)
}

#' Function imported from https://github.com/th1vairam/CovariateAnalysis
#' Modified from http://stackoverflow.com/questions/13088770/
#' Function to find linearly dependednt columns of a matrix
#' @param mat an input design matrix.
#' @return a list of independent columns
#' @examples
#' data("iris")
#' designMat <- getDesignMatrix(iris)
#' linColumnFinder(designMat$design)
#' @export
linColumnFinder <- function(mat) {
    mat[is.na(mat)] <- 0

    # If the matrix is full rank then we're done
    if (qr(mat)$rank == ncol(mat)) {
        return(list(
            indepCols = seq(1, ncol(mat), 1),
            relations = "Matrix is of full rank"
        ))
    }

    m <- ncol(mat)
    # cols keeps track of which columns are linearly independent
    cols <- 1
    allMessages <- c()
    for (i in seq(2, m)) {
        ids <- c(cols, i)
        mymat <- mat[, ids]
        if (qr(mymat)$rank != length(ids)) {
            # Regression the column of interest on the previous columns to
            # figure out the relationship
            o <- stats::lm(mat[, i] ~ as.matrix(mat[, cols]) + 0)
            # Construct the output message
            start <- paste0(colnames(mat)[i], " = ")
            # Which coefs are nonzero
            nz <- !(abs(stats::coef(o)) <= .Machine$double.eps^0.5)
            tmp <- colnames(mat)[cols[nz]]
            vals <- paste(stats::coef(o)[nz], tmp, sep = "*", collapse = " + ")
            amessage <- paste0(start, vals)
            allMessages <- c(allMessages, amessage)
        } else {
            # If the matrix subset was of full rank
            # then the newest column in linearly independent
            # so add it to the cols list
            cols <- ids
        }
    }
    return(list(indepCols = cols, relations = allMessages))
}

####
#
# Helper functions for differentialPathwayAnalysis
#
####



#' function to get residuals with respect to a set of covariates
#'
#' @param covariates a covariate dataframe.
#' @param adjustCovars covariates to adjust for
#' @param pathSumStats an expression matrix
#' @return a matrix of adjusted expression
getResidual <- function(covariates,
                        adjustCovars,
                        pathSumStats) {
    # Getting pathway residuals
    resDesingMat <- getDesignMatrix(covariates[, c(adjustCovars), drop = FALSE],
        intercept = FALSE)

    resDesingMat$design <- resDesingMat$design[,
        linColumnFinder(resDesingMat$design)$indepCols]

    fitResiduals <- limma::lmFit(pathSumStats, resDesingMat$design)
    fitResiduals <- limma::residuals.MArrayLM(fitResiduals, pathSumStats)

    return(fitResiduals)
}

#' function to get a DE table
#'
#' @param expMat an expression matrix
#' @param designMat a design Matrix
#' @param contrastsName the contrast to perform
#' @return a table of differential expression
getDiffExpTable <- function(expMat,
                            designMat,
                            contrastsName) {
    fits <- limma::lmFit(expMat, designMat$design)

    colnames(fits$coefficients) <- gsub("-", "_", colnames(fits$coefficients))

    contrast <- limma::makeContrasts(contrasts = c(contrastsName),
                                    levels = colnames(fits$coefficients))

    fitContrast <- limma::contrasts.fit(fits, contrasts = contrast)
    fitContrast <- limma::eBayes(fitContrast)
    tT <- limma::topTable(fitContrast, adjust = "fdr", sort.by = "p",
                    number = Inf)
    tT$contrast <- contrastsName

    return(tT)
}


.constrastHelper <- function(covariates, condition, contrastConds) {
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
    return(contrastsName)
}

.pathwayCleaner <- function(pathways, id, geneCounts, minPathSize) {
    # select pathways with genes in the count data and a minimum pathway size
    pathways <- as.data.frame(pathways)
    genesPathways <- pathways[pathways[, id] %in% rownames(geneCounts), ]

    genesPathways <- genesPathways |> dplyr::group_by(Pathway) |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::filter(n >= minPathSize)

    pathways <- pathways |>
        dplyr::filter(Pathway %in% genesPathways$Pathway)
    return(pathways)
}

.qNormalizer <- function(quantileNorm, pathSumStats, covariates) {

    if (quantileNorm == TRUE) {
        pathwayNames <- rownames(pathSumStats)
        pathSumStats <-
            preprocessCore::normalize.quantiles(pathSumStats)

        rownames(pathSumStats) <- pathwayNames
        colnames(pathSumStats) <- rownames(covariates)
    }
    return(pathSumStats)
}


.dirChecks <- function(outDir) {
    if (substring(outDir, nchar(outDir)) != "/") {
        outDir <- paste0(outDir, "/")
    }
    if (!dir.exists(outDir)) {
        stop("Output directory does not exist.")
    }
    return(outDir)
}

.pathSumCleaner <- function(pathSumStats) {
    # filter pathways with na values in the pathway summary statistics
    pathSumStats <- pathSumStats[rowSums(is.na(pathSumStats)) == 0, ]
    pathSumStats <- apply(pathSumStats, 2, function(x) {
        (x - mean(x)) / stats::sd(x)
    })
    return(pathSumStats)
}

.tScoreMaker <- function(deGenes, id) {
    tScores <- NULL
    if (!is.null(deGenes)) {
        tScores <- deGenes |>
            dplyr::mutate(!!id := rownames(deGenes)) |>
            dplyr::select(c(ENSEMBL, t))
    }
    return(tScores)
}
####
#
# Helper functions for miRNAPathwayEnrichment
#
####



#' function to align a list of sets and a reference universe
#'
#' @param pathwaySets a list of sets
#' @param universe all set  elements must be a subset of universe
#' @return a list of sets, aligned to universe
alignToUniverse <- function(pathwaySets,
                            universe) {
    pathwaySets <- lapply(pathwaySets, function(x) {
        x[x %in% universe]
    })
    return(pathwaySets)
}


#' Pairwise enrichment analysis between two given lists of sets
#'
#' @param mirSets a list of targets of miRNAs
#' @param pathwaySets a list of pathways
#' @param pathsRef universe of genes.
#' @param numCores number of cores to calculate the results.
#' @return enrichment analysis results
enrichAllPairs <- function(mirSets,
                            pathwaySets,
                            pathsRef,
                            numCores) {
    iterator <- (merge(names(mirSets), names(pathwaySets)))
    iterator <- iterator |> dplyr::mutate_all(as.character)
    all <- length(pathsRef)

    # find enrichment p-value of each miRNA target set and each pathway set
    enrichs <- parallel::mclapply(seq_len(nrow(iterator)),
                    function(y) {
                        x <- iterator[y, ]
                        q <- length(intersect(unlist(pathwaySets[x[[2]]]),
                                            unlist(mirSets[x[[1]]])))
                        m <- length(unlist(mirSets[x[[1]]]))
                        n <- all - m
                        k <- length(unlist(pathwaySets[x[[2]]]))
                        pval <- stats::phyper(
                            q - 1, m, n, k, lower.tail = FALSE, log.p = FALSE)

                        return(c("pval" = pval, "Intersect" = q,
                            "mirset_Size" = m,
                            "not_mirset" = n,
                            "pathway_Size" = k,
                            "ratio_in" = q / (k - q),
                            "ratio_out" = (m - q) / n,
                            "ratio_ratios" = (q / (k - q)) / ((m - q) / n)))
                        }, mc.cores = numCores)

    temp <- do.call(rbind, enrichs)
    temp <- as.data.frame(temp)
    iterator <- cbind(iterator, temp)

    return(iterator)
}

####
#
# Helper functions for mappingPathwayClusters
#
####


#' Creates a network out of pcxn table
#'
#' @param pcxn pathways network edge list of pathways
#' @param edgeFDR FDR threshold for pathway-pathway adjusted p-values;
#'   filter edges with adjusted p-values less than given threshold
#' @param correlationCutOff cut-off threshold for pathway-pathway correlation;
#'   filter pathways with correlation less than given threshold
#' @param weighted True if you wish to include correlation weights in clustering
#' @return enrichment analysis results
pcxnToNet <- function(pcxn,
                        edgeFDR,
                        correlationCutOff,
                        weighted) {
    pcxn$pAdjust <- stats::p.adjust(pcxn$p.value, method = "fdr")
    res <- pcxn[pcxn$pAdjust < edgeFDR, ]
    res <- res[abs(res$PathCor) > correlationCutOff, ]
    allPathways <- union(unique(res$Pathway.B), unique(res$Pathway.A))
    net <- res[, c("Pathway.A", "Pathway.B", "PathCor")]
    if (weighted) {
        net$weight <- abs(net$PathCor)
    }
    net <- as.matrix(net)
    net <- igraph::graph_from_data_frame(net, directed = FALSE)

    igraph::set_vertex_attr(net, "col", value = "red")
    return(list("net" = net, "allPathways" = allPathways))
}

.levelFixerMapper <- function(clusts) {
    zz <- as.factor(clusts$membership)
    zz <- forcats::fct_infreq(zz, ordered = NA)
    levels(zz) <- seq_along(unique(zz))
    clusts$membership <- as.numeric(zz)
    return(clusts)
}

.shapeColNet <- function(subNet, dePathways1) {
    igraph::E(subNet)$color <-
        ifelse(as.numeric(igraph::E(subNet)$PathCor) > 0, "#E41A1C", "#377EB8")

    shapeIndex <- rownames(dePathways1) %in% igraph::V(subNet)$name
    shapeIndex <- dePathways1[shapeIndex, c("logFC", "AveExpr")]
    shapeIndex <- shapeIndex[igraph::V(subNet)$name, ]
    shapeIndex <- shapeIndex$logFC > 0
    igraph::V(subNet)$shape   <- ifelse(shapeIndex, "square", "circle")
    return(subNet)
}


####
#
# Helper functions for prioritizeMicroRNAs
#
####

.checkAddressDirs <- function(outDir, saveCSV, dataDir, saveSampling) {
    if (!dir.exists(outDir) && saveCSV) {
        warning("Output directory does not exist.")
        dir.create(outDir, recursive = TRUE)
    }
    if (!dir.exists(dataDir) && saveSampling) {
        stop("Data directory does not exist.")
    }
}

.cleanEnrichInput <- function(enriches0, enrichmentFDR) {
    enriches <- enriches0 |> dplyr::filter(Intersect != 0)
    enriches <- enriches |> dplyr::group_by(y) |>
        dplyr::mutate(path_fdr = stats::p.adjust(pval, method = "fdr"))
    enriches <- enriches |>
        dplyr::mutate(hit = ifelse(path_fdr < enrichmentFDR, 1, 0))
    return(enriches)
}


.makeSelector <- function(enriches, pathways) {
    # Count miRNA-pathway enrichment
    tempEnrich <- enriches[enriches$y %in% pathways, ]
    selector <- tempEnrich |>
        dplyr::group_by(x) |>
        dplyr::summarise("cluster_hits" = sum(hit))
    return(selector)
}


.saveCsvWriter <- function(saveCSV, prefix, sampRate, clustNo, selector,
                        outDir) {
    if (saveCSV == TRUE) {
        saveName <- paste0(prefix, sampRate, "_samples_clustNo_", clustNo,
                            ".csv")
        sayThis <- paste0(saveName, " saved!")
        message(sayThis)
        utils::write.csv(selector, paste0(outDir, saveName))
    }
}


.mkSampleDatDir <- function(dataDir, prefix, saveSampling, m, sampRate) {
    samplingDataDir <- paste0(dataDir, prefix, "Sampling_Data/")

    if ((saveSampling == TRUE) && (!dir.exists(samplingDataDir))) {
        dir.create(samplingDataDir, recursive = TRUE)
    }
    tempName <- paste0(prefix, m, "_", sampRate, "_samples.RDS")
    samplingDataFile <- paste0(samplingDataDir, "/", tempName)

    return(samplingDataFile)
}

.mSelectorMaker <- function(m, enrichNull, mSelector, sampRate, fn, nPaths,
                        sampDatFile, saveSampling, numCores, autoSeed, coverFn,
                        runJackKnife, pathways) {

    if (m %in% c("aggInv", "aggLog")) {
        samplingData <- samplingDataBase(enrichNull, mSelector,
                        sampRate, fn, nPaths, sampDatFile, jackKnife = FALSE,
                        saveSampling, numCores = numCores, autoSeed = autoSeed)
        mSelector <- methodProbBase(
            samplingData[[paste0("SampSize_", 100)]], mSelector,
            m = m, nPaths = nPaths, coverFn = coverFn)
    } else {
        names(mSelector)[2:3] <- paste0(m, "_", names(mSelector)[2:3])
    }
    sayThis <- paste0(m, " Method Done")
    message(sayThis)
    # perform jack-knife

    if ((runJackKnife == TRUE) && (m != "sumz") && (m != "sumlog")) {

        samplingData <- samplingDataBase(enrichNull, mSelector, sampRate, fn,
                    nPaths, sampDatFile, jackKnife = FALSE, saveSampling,
                    numCores = numCores, autoSeed = autoSeed)

        mSelector <- jackKnifeBase(selector = mSelector, pathways = pathways,
                    enrichNull = enrichNull, fn = fn,
                    jackKnifeData =  samplingData[[paste0("SampSize_", 100)]],
                    m = m, numCores = numCores)

        sayThis <- paste0(m, " JackKnifing Method Done!")
        message(sayThis)

        mSelector <- mSelector[, c(1, 3, 2, 4)]
    } else {
        mSelector <- mSelector[, c(1, 3, 2)]
    }
    return(mSelector)
}


.expandSelector <- function(selector, mSelector, m) {

    selector <- merge(selector, mSelector, all = TRUE)
    met <- paste0(m, "_pval")
    selector <- selector |> dplyr::arrange(!!rlang::sym(met))
    met2 <- paste0(m, "_fdr")
    selector <- selector |> dplyr::mutate(
        !!met2 := stats::p.adjust(p = !!rlang::sym(met), method = "fdr")
    )
    return(selector)
}
####
#
# Helper functions for clusterPlotter
#
####

.clusterPlotHelper <- function(plotSave, figDir, subNet, nodeColors,
                                legendCats) {
    if (plotSave) {
        grDevices::pdf(paste0(figDir, "PCxNCorGraph.pdf"),
                    width = 18, height = 11)
    }
    plot(subNet, vertex.size = 5, vertex.label = NA, vertex.color = nodeColors)
    graphics::legend(x = "bottomleft", legend = legendCats$attr, pch = c(0, 1),
                    bty = "n", cex = 1.6)
    graphics::legend(x = "topleft", legend = c("Positive Cor", "Negative Cor"),
                    col = c("#E41A1C", "#377EB8"), lty = 1, lwd = 2,
                    cex = 1.6, bty = "n")
    if (plotSave)
        grDevices::dev.off()
}

.clusterSubPlotHelper <- function(subplot, topClusters, clstMems, subNet,
                                figDir, cols, legendCats) {
    if (subplot == TRUE) {
        for (k in seq_len(topClusters)) {
            keep <- which((clstMems) == k)
            subNet2 <- igraph::induced_subgraph(subNet, keep)
            if (length(igraph::V(subNet2)) < 2) next

            grDevices::pdf(paste0(figDir, "PCxNCorGraph_",
                            "Cluster_", k, ".pdf"))
            plot(subNet2, edge.width = 1.3, vertex.size = 5, vertex.label = NA,
                vertex.color = cols[clstMems[keep]], legend = TRUE,
                layout = igraph::layout.fruchterman.reingold)
            graphics::legend(x = "bottomleft", legend = legendCats$attr,
                            pch = c(0, 1), bty = "n", cex = 1.4)
            graphics::legend(x = "topleft",
                            legend = c("Positive Cor", "Negative Cor"),
                            col = c("#E41A1C", "#377EB8"),
                            lty = 1, lwd = 2, cex = 1.4, bty = "n")
            grDevices::dev.off()
        }
    }
}
