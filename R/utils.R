utils::globalVariables(c(
    "%<>%", ".", ":=", "ENSEMBL", "ES", "ES2", "Intersect", "Pathway", "hit",
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
#' @param outDir Address to save an RDS for a table of pathway-gene association
#' @return pathExpTab Table of pathway-gene association.
#' @export
pathwayGeneTab <- function(pathAdress = NA,
                            pathwayList = NA,
                            outDir = NA) {
    # Development notes: make compatible with direct data input
    #       check and throw errors if the address is valid
    #       make it work for all anotation types of genes Eg. Symbol or ENSEMBL.
    #       Make a  pathwayGeneTab going out with the package.

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
            clusterProfiler::bitr(
                x,
                "ENTREZID",
                "ENSEMBL",
                OrgDb = org.Hs.eg.db::org.Hs.eg.db
            )
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
#' @export

pathwaySummary <- function(exprsMat,
                            pathwayRef,
                            id = "ENSEMBL",
                            zNormalize = FALSE,
                            method = FALSE,
                            deGenes = NULL,
                            trim = 0,
                            tScores = NULL) {

    # There is a confusion about the data format here.
    # Make sure it is consistent.
    # Current version only works with ENSEMBL.
    exprsMat <- tibble::rownames_to_column(as.data.frame(exprsMat), var = id)

    if (!is.null(deGenes)) {
        if (is.null(tScores)) {
            stop("Provide tscores/pvalues")
        }

        pathwayRef <- dplyr::inner_join(pathwayRef, tScores, by = c("ENSEMBL")
        )

        pathwayRef %<>% dplyr::group_by(., Pathway) %>%
            dplyr::filter(., abs(t) >= stats::median(abs(t))) %>%
            dplyr::select(., -t)
    }

    if (method == "none") {
        exprsMat <- exprsMat %>%
            dplyr::mutate_if(., is.numeric, function(x) {
                x
                })
    } else if (method == "x") {
        exprsMat <- exprsMat %>%
            dplyr::mutate_if(., is.numeric, rank) %>%
            dplyr::mutate_if(., is.numeric, function(x) {
                x
                })
    } else if (method == "x2") {
        exprsMat <- exprsMat %>%
            dplyr::mutate_if(., is.numeric, rank) %>%
            dplyr::mutate_if(., is.numeric, function(x) {
                x * x
            })
    } else {
        stop("invalid choice of summarization function")
    }

    pathExpTab <- dplyr::inner_join(pathwayRef, exprsMat, by = id)

    pathExpTab <- pathExpTab %>%
        dplyr::group_by(., Pathway) %>%
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
pCutFn <- function(enriches, pathways, isSelector, thresh = 0.05) {
    if (isSelector == TRUE) {
        enriches <- enriches %>%
            dplyr::mutate(., hit2 = ifelse(pval < thresh, 1, 0))
    }

    tempEnrich <- enriches[enriches$y %in% pathways, ]
    selector <- tempEnrich %>%
        dplyr::group_by(x) %>%
        dplyr::summarise(n = n(), k = sum(hit2)) %>%
        dplyr::arrange(., x)

    if (isSelector == TRUE) {
        selector <- selector %>% dplyr::filter(., k > thresh * length(pathways))
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
aggInvFn <- function(enriches, pathways, isSelector = TRUE, thresh = NULL) {
    if (isSelector == TRUE) {
        enriches <- enriches %>% dplyr::mutate(., ES2 = stats::qnorm(1 - pval))
        minES <- min(enriches$ES2[!is.infinite(enriches$ES2)])
        enriches <- enriches %>%
            dplyr::mutate(., ES2 = ifelse(is.infinite(.$ES2), minES, .$ES2))
    }

    tempEnrich <- enriches[enriches$y %in% pathways, ]
    selector <- tempEnrich %>%
        dplyr::group_by(x) %>%
        dplyr::summarise(n = n(), k = mean(ES2)) %>%
        dplyr::arrange(., x)

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
aggLogFn <- function(enriches, pathways, isSelector, thresh = 0) {
    enriches <- enriches %>% dplyr::mutate(., ES = -log(pval))
    if (isSelector == TRUE) {
        enriches <- enriches %>% dplyr::mutate(., ES = -log(pval))
    }

    tempEnrich <- enriches[enriches$y %in% pathways, ]
    selector <- tempEnrich %>%
        dplyr::group_by(x) %>%
        dplyr::summarise(n = n(), k = mean(ES))

    if (isSelector == TRUE) {
        selector <- selector %>%
            dplyr::filter(., k * n >= thresh * length(pathways))
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
sumzFn <- function(enriches, pathways, isSelector, thresh = NULL) {
    enriches1 <- enriches %>%
        dplyr::mutate(., pval = ifelse(pval >= 0.999, 0.999, pval))

    enriches1 <- enriches1 %>%
        dplyr::mutate(., pval = ifelse(pval <= 1.0e-16, 1.0e-16, pval))

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
sumlogFn <- function(enriches, pathways, isSelector, thresh = NULL) {
    enriches1 <- enriches %>% dplyr::mutate(., pval = ifelse(pval >= 0.999,
                                    0.999, pval))
    enriches1 <- enriches1 %>% dplyr::mutate(., pval = ifelse(pval <= 1.0e-16,
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
pCutCoverFn <- function(selector, coverName) {
    selector <- selector %>%
        dplyr::mutate(., !!coverName := k / n)
    return(selector)
}


#' Internal function for modification of prioritization.
#' @param selector a prioritzation table
#' @param coverName a new column name
#' @return an updated scoring of miRNAs in a cluster of pathways
aggInvCoverFn <- function(selector, coverName) {
    selector <- selector %>%
        dplyr::mutate(., !!coverName := k)
    return(selector)
}

#' Internal function for modification of prioritization.
#' @param selector a prioritzation table
#' @param coverName a new column name
#' @return an updated scoring of miRNAs in a cluster of pathways
aggLogCoverFn <- aggInvCoverFn

#' Internal function for modification of prioritization.
#' @param selector a prioritzation table
#' @param coverName a new column name
#' @return an updated scoring of miRNAs in a cluster of pathways
sumzCoverFn <- aggInvCoverFn

#' Internal function for modification of prioritization.
#' @param selector a prioritzation table
#' @param coverName a new column name
#' @return an updated scoring of miRNAs in a cluster of pathways
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
samplingDataBase <- function(enrichNull,
                            selector,
                            sampRate,
                            fn,
                            nPaths,
                            samplingDataFile,
                            jackKnife = FALSE,
                            saveSampling,
                            numCores = 1,
                            autoSeed = TRUE) {
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

                selNull <- fn(
                    enriches = enrichNull,
                    pathways = nullPaths,
                    isSelector = FALSE
                )
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
            print(paste0(samplingDataFile, " saved."))
        }
    } else {
        print(paste0("Skipping sampling, ", samplingDataFile, " exists."))
        outList <- readRDS(samplingDataFile)
        print(paste0(samplingDataFile, " loaded."))
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
    selector <- selector %>%
        dplyr::mutate(., !!pvalName := pVals) %>%
        coverFn(., coverName) %>%
        dplyr::select(., -c(k, n))
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
    selector <- selector %>%
        dplyr::mutate(., !!jackKnifeName := means)

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
    rowNamesTemp <- rownames(covariatesDataFrame)
    colNamesTemp <- colnames(covariatesDataFrame)

    factorCovariateNames <-
        names(covariatesDataFrame)[vapply(covariatesDataFrame,
                                        is.factor,
                                        logical(1))]

    factorCovariateNames <-
        setdiff(
            factorCovariateNames,
            factorCovariateNames[
                !(factorCovariateNames %in% colnames(covariatesDataFrame))
            ]
        )

    numericCovariateNames <- setdiff(colNamesTemp, factorCovariateNames)

    # Ensure the factors are in fact of type factor, and the quantitative
    # variables are numeric:
    covariatesDataFrame <-
        as.data.frame(lapply(colnames(covariatesDataFrame), function(column) {
            if (column %in% factorCovariateNames) {
                fac <- as.factor(covariatesDataFrame[, column])
                if (column %in% names(reLevels)) {
                    fac <- stats::relevel(fac, ref = reLevels[[column]])
                }
                return(fac)
            } else {
                return(as.numeric(covariatesDataFrame[, column]))
            }
        }))

    rownames(covariatesDataFrame) <- rowNamesTemp
    colnames(covariatesDataFrame) <- colNamesTemp

    contra <- NULL
    maxNumCat <- Inf
    catData <- covariatesDataFrame[, factorCovariateNames, drop = FALSE]
    if (ncol(catData) > 0) {
        numCats <- vapply(
            colnames(catData),
            function(col) nlevels(factor(catData[, col])),
            numeric(1)
        )

        excludeCategoricalCols <-
            names(numCats)[numCats <= 1 | numCats > maxNumCat]

        if (!is.null(excludeCategoricalCols) &&
            length(excludeCategoricalCols) > 0) {
            warning(paste("Excluding categorical variables with less than 2",
                ifelse(is.infinite(maxNumCat), "",
                        paste(" or more than ", maxNumCat, sep = "")),
                " categories: ", paste(paste("'", excludeCategoricalCols,
                                    "'", sep = ""), collapse = ", "),
                sep = ""))

            factorCovariateNames <- setdiff(
                factorCovariateNames,
                excludeCategoricalCols
            )
            covariatesDataFrame <-
                covariatesDataFrame[,
                !(colnames(covariatesDataFrame) %in% excludeCategoricalCols),
                drop = FALSE]
        }

        # Inspired by http://stackoverflow.com/questions/4560459/
        #
        # And, already ensured above that
        # covariatesDataFrame[, factorCovariateNames] satisfies:
        # 1) fac is of type factor.
        # 2) fac is releveled as designated in reLevels.
        if (intercept) {
            contra <- lapply(
                factorCovariateNames,
                function(column) {
                    fac <- covariatesDataFrame[, column]
                    fac <- stats::contrasts(fac)
                }
            )
        } else {
            contra <- lapply(
                factorCovariateNames,
                function(column) {
                    fac <- covariatesDataFrame[, column]
                    fac <- stats::contrasts(fac, contrasts = FALSE)
                }
            )
        }
        names(contra) <- factorCovariateNames
    }

    # Inspired by http://stackoverflow.com/questions/5616210/
    currentNAAction <- getOption("na.action")
    # Model matrix will now include "NA":
    options(na.action = "na.pass")

    if (intercept) {
        design <- stats::model.matrix(~., data = covariatesDataFrame,
                            contrasts.arg = contra
        )
    } else {
        design <- stats::model.matrix(~ 0 + ., data = covariatesDataFrame,
                            contrasts.arg = contra
        )
    }

    rownames(design) <- rownames(covariatesDataFrame)

    options(na.action = currentNAAction)

    return(list(
        design = design,
        covariates = colNamesTemp,
        factorsLevels = lapply(contra, colnames),
        numericCovars = numericCovariateNames,
        covariatesDataFrame = covariatesDataFrame
    ))
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
            message <- paste0(start, vals)
            allMessages <- c(allMessages, message)
        } else {
            # If the matrix subset was of full rank
            # then the newest column in linearly independent
            # so add it to the cols list
            cols <- ids
        }
    }
    return(list(indepCols = cols, relations = allMessages))
}
