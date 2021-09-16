# 
# sampling.data.base <- function(enrich.null,
#                                selector,
#                                samp.rate,
#                                fn,
#                                n_paths,
#                                sampling.data.file,
#                                jack.knife=FALSE,
#                                save.sampling,
#                                num.cores=1){
#   if(!all(hasName(selector,c("x"))))
#     stop("The selector table needs a column x (miRNA name)")
#   if(!all(hasName(enrich.null,c("x","y","pval"))))
#     stop(paste0("The enrichment table needs a column x (miRNA name)",
#                 ",a column y (pathway name), and a pval column"))
#   
#   if (!file.exists(sampling.data.file)){
#     all.paths   <- unique(enrich.null$y)
#     temp.n_paths <- n_paths
#     
#     if (jack.knife==TRUE) {temp.n_paths <- n_paths-1}
#     
#     # build null distribution of K 
#     temp   <- mclapply(1:(samp.rate), function(Y){
#       set.seed(Y)        
#       null.paths <- sample(all.paths,temp.n_paths,replace = F)
#       sel.null   <- fn(enriches=enrich.null, pathways=null.paths, is.selector=F)
#       return(sel.null$k)
#       
#     }, mc.cores = num.cores )
#     
#     
#     temp <- do.call(rbind, temp)
#     temp <-  t(temp)
#     
#     rownames(temp) <- selector$x
#     colnames(temp) <- sapply(1:(samp.rate), function(Y){paste0("sample_", Y)})
#     
#     if (save.sampling==TRUE) {
#       saveRDS(temp, file=sampling.data.file)
#       print(paste0(sampling.data.file, " saved."))
#     }
#     
#   }
#   else{
#     print(paste0("Skipping sampling, ", sampling.data.file, " exists."))
#     temp <- readRDS(sampling.data.file)
#     print(paste0(sampling.data.file, " loaded."))
#   }
#   
#   return(temp)
#   
# }




# method.prob.base <- function(sampling.data, selector, m, cover.fn=NULL){
#   
#   if(!all(hasName(selector,c("x","k"))))
#     stop("The selector table needs a column x (miRNA name) and a column k (miRNA hits)")
#   
#   # obtain means and sds for distribution, assume CLT
#   means <- rowMeans(sampling.data) 
#   sds <- apply(sampling.data, 1, sd)
#   
#   pval.name <- paste0(m, '_pval')
#   cover.name <- paste0(m, '_cover')
#   
#   # obtain p-vals
#   p_vals <- pnorm(selector$k, mean=means, sd=sds, lower.tail=FALSE)
#   selector <- selector %>%
#     dplyr::mutate(., !!pval.name := p_vals) %>%
#     cover.fn(., cover.name) %>%
#     dplyr::select(.,-c(k,n))
#   return(selector)
# }
# 
# 


method.prob.base3 <- function(sampling.data, selector,name.tags = "", method.tags="",
                              n_paths = 100, cover.fn=NULL){
    
    if(!all(hasName(selector,c("x","k"))))
        stop("The selector table needs a column x (miRNA name) and a column k (miRNA hits)")
    
    # obtain means and sds for distribution, assume CLT
    means      <- rowMeans(sampling.data) 
    sds        <- apply(sampling.data, 1, sd)
    sds        <- sds *10/sqrt(n_paths)
    pval.name  <- paste0('pval', method.tags)
    cover.name <- paste0(method.tags, '_cover')
    
    # obtain p-vals
    p_vals <- pnorm(selector$k, mean=means, sd=sds, lower.tail=FALSE)
    selector <- selector %>%
        dplyr::mutate(., !!pval.name := p_vals) %>%
        cover.fn(., cover.name) %>%
        dplyr::select(.,-c(k,n))
    
    colnames(selector)[-1] <- paste0(name.tags,colnames(selector)[-1])
    
    return(selector)
}






jack.knife.base <- function(selector, pathways, enrich.null, fn, jack.knife.data, m, num.cores=1){
    
    # obtain means and sds for distribution, assume CLT
    sample.means <- rowMeans(jack.knife.data) 
    sample.sds <- apply(jack.knife.data, 1, sd)
    
    # remove one pathway at a time and obtain K for each miRNA 
    temp1 <- mclapply(1:length(pathways), function(X){
        temp.pathways  <- pathways[-X]
        temp.selector   <- fn(enriches=enrich.null, pathways=temp.pathways, is.selector=F)
        # obtain p-values using the means and sds obtain above
        p_vals <- pnorm(temp.selector$k, mean=sample.means, sd=sample.sds, lower.tail = FALSE)
        return(p_vals)
    }, mc.cores = num.cores)
    
    
    # rows <- number of pathways; col <- number of miRNAs
    # contain p-values
    temp1 <- do.call(rbind, temp1)
    temp1 <- t(temp1)
    
    ## added from here
    
    jack.knife.name <- paste0(m, "_pval_jk")
    
    # obtain aggregate p-values
    means <- rowMeans(temp1)
    selector <- selector %>%
        dplyr::mutate(., !!jack.knife.name := means)
    
    return(selector)
}



miRNAPrioritization <- function(enriches0,
                                pathway.clusters,
                                method='AggInv',
                                method.thresh=NULL,
                                mir.path.fdr.thresh=0.25,
                                top.clust=2,
                                samp.rate=1000,
                                out.dir='',
                                data.dir='',
                                save.sampling=T,
                                run.jack.knife=T,
                                save.jack.knife=F,
                                num.cores=1,
                                save.csv=T,
                                prefix='')
{
    
    if (substring(out.dir, nchar(out.dir))!='/')
        out.dir <- paste0(out.dir, '/')
    if (!dir.exists(out.dir))
        stop('Output directory does not exist.')
    
    if(out.dir=='/')
        out.dir=''
    
    if (substring(data.dir, nchar(data.dir))!='/')
        data.dir <- paste0(data.dir, '/')
    if (!dir.exists(data.dir))
        stop('Data directory does not exist.')
    
    if(data.dir=='/')
        data.dir=''
    
    output = list()
    
    # count miRNA-pathway enrichment with p-value less than threshold
    enriches   <- enriches0 %>% filter(., Intersect != 0)
    enriches   %<>% group_by(., y) %>% mutate(., path_fdr = p.adjust(pval, method = "fdr"))
    enriches   <- enriches %>% mutate(.,hit =ifelse(path_fdr < mir.path.fdr.thresh, 1, 0))
    
    for (clustNo in 1:top.clust){
        
        clustName <- paste0('Cluster', clustNo)
        
        print(paste0('Working on ', clustName, '.'))
        
        # select pathways in cluster
        pathways    <- as.character(pathway.clusters[pathway.clusters$cluster == clustNo, ]$Pathway)
        n_paths     <- length(pathways)
        
        # formulate number of miRNA-pathway enrichment with p-value less than threshold for each miRNA
        temp.enrich <- enriches[enriches$y %in% pathways, ]
        selector   <- temp.enrich %>% group_by(x)  %>%
            dplyr::summarise(., "cluster_hits" = sum(hit))
        
        # perform p-value aggregation based on methodlogy provided
        for (i in 1:length(method)){
            m <- method[i]
            
            print(paste0('Performing ', m, ' function.'))
            
            fn       <- get(paste0(m, '.fn'))
            cover.fn <- get(paste0(m, '.cover.fn'))
            
            if (!is.null(method.thresh)){
                m.thresh <- method.thresh[i]
                temp <- fn(enriches=enriches0, pathways, is.selector=T, thresh=m.thresh)
            } else{
                temp <- fn(enriches=enriches0, pathways, is.selector = T)
            }
            
            m.selector <- temp$selector
            m.enriches0 <- temp$enriches0
            
            if (nrow(m.selector) < 3){
                print(paste0('Skipping ', m, ' function due to low number of miRNA after filter'))
                next
            }
            
            enrich.null <- m.enriches0 %>% filter(.,x %in% m.selector$x)
            
            sampling.data.dir <- paste0(data.dir, prefix, 'Sampling_Data/', clustName)
            
            if (save.sampling==TRUE){
                if(!dir.exists(sampling.data.dir))dir.create(sampling.data.dir,recursive = T)
            }
            
            sampling.data.filename  <- paste0(prefix, m, "_", samp.rate, "_samples.RDS")
            sampling.data.file <- paste0(sampling.data.dir, '/', sampling.data.filename)
            
            # perform sampling
            sampling.data <- sampling.data.base(enrich.null, m.selector, samp.rate, fn, n_paths, sampling.data.file, jack.knife=FALSE, save.sampling = save.sampling)
            
            m.selector <- method.prob.base(sampling.data, m.selector, m, cover.fn)
            
            print(paste0(m, " Method Done"))
            
            # perform jack-knife
            if (run.jack.knife ==T){
                jack.knife.data.filename  <- paste0(prefix, m, "_", samp.rate, "_samples_jack_knifed.RDS")
                jack.knife.data.file <- paste0(sampling.data.dir, '/', jack.knife.data.filename)
                
                jack.knife.data <- sampling.data.base(enrich.null, m.selector, samp.rate, fn, n_paths, 
                                                      sampling.data.file=jack.knife.data.file, jack.knife=TRUE, save.sampling = save.jack.knife)
                
                m.selector <- jack.knife.base(m.selector, pathways, enrich.null, fn, jack.knife.data, m)
                
                print(paste0(m, " JackKnifing Method Done!"))
                
                m.selector <- m.selector[,c(1,3,2,4)]
            } else {
                m.selector <- m.selector[,c(1,3,2)]
            }
            
            selector <- merge(selector, m.selector, all=T)
        }
        
        if (save.csv==T){
            save.name <- paste0(prefix, samp.rate, '_samples_clustNo_', clustNo, '.csv')
            print(paste0(save.name, ' saved!'))
            write.csv(selector, paste0(out.dir, save.name))
        }
        
        output[[clustName]] <- selector
    }
    
    return(output)
}






library(TCGAbiolinks)
library(CovariateAnalysis)
library(reshape2)
library(ggplot2)
library(psych)

TCGAbatch_Correction_revised <- function (tabDF,
                                          batch.factor = NULL,
                                          adjustment = NULL,
                                          ClinicalDF = data.frame(),
                                          UnpublishedData = FALSE,
                                          AnnotationDF = data.frame(),
                                          plot = TRUE)
    
{
    if( UnpublishedData == TRUE) {
        batch.factor <- as.factor(AnnotationDF$Batch)
        batch_corr <- sva::ComBat(dat = tabDF, batch = batch.factor, par.prior = TRUE, prior.plots = plot)
    }
    
    if( UnpublishedData == FALSE) {
        
        if (length(batch.factor) == 0 & length(adjustment) == 0)
            message("batch correction will be skipped")
        else if (batch.factor %in% adjustment) {
            stop(paste0("Cannot adjust and correct for the same factor|"))
        }
        my_IDs <- get_IDs(tabDF)
        if (length(batch.factor) > 0 || length(adjustment) > 0)
            if ((nrow(ClinicalDF) > 0 & batch.factor == "Year") ||
                ("Year" %in% adjustment == TRUE & nrow(ClinicalDF) >
                 0)) {
                names(ClinicalDF)[names(ClinicalDF) == "bcr_patient_barcode"] <- "patient"
                ClinicalDF$age_at_diag_year <- floor(ClinicalDF$age_at_diagnosis/365)
                ClinicalDF$diag_year <- ClinicalDF$age_at_diag_year +
                    ClinicalDF$year_of_birth
                diag_yearDF <- ClinicalDF[, c("patient", "diag_year")]
                Year <- merge(my_IDs, diag_yearDF, by = "patient")
                Year <- Year$diag_year
                Year <- as.factor(Year)
            }
        else if (nrow(ClinicalDF) == 0 & batch.factor == "Year") {
            stop("Cannot extract Year data. Clinical data was not provided")
        }
        Plate <- as.factor(my_IDs$plate)
        Condition <- as.factor(my_IDs$condition)
        TSS <- as.factor(my_IDs$tss)
        Portion <- as.factor(my_IDs$portion)
        Sequencing.Center <- as.factor(my_IDs$center)
        design.matrix <- model.matrix(~Condition)
        design.mod.combat <- model.matrix(~Condition)
        options <- c("Plate", "TSS", "Year", "Portion", "Sequencing Center")
        
        if (length(batch.factor) > 1)
            stop("Combat can only correct for one batch variable. Provide one batch factor")
        if (batch.factor %in% options == FALSE)
            stop(paste0(o, " is not a valid batch correction factor"))
        
        for (o in adjustment) {
            if (o %in% options == FALSE)
                stop(paste0(o, " is not a valid adjustment factor"))
            
        }
        adjustment.data <- c()
        for (a in adjustment) {
            if (a == "Sequencing Center")
                a <- Sequencing.Center
            adjustment.data <- cbind(eval(parse(text = a)), adjustment.data)
        }
        if (batch.factor == "Sequencing Center")
            batch.factor <- Sequencing.Center
        batchCombat <- eval(parse(text = batch.factor))
        if (length(adjustment) > 0) {
            adjustment.formula <- paste(adjustment, collapse = "+")
            adjustment.formula <- paste0("+", adjustment.formula)
            adjustment.formula <- paste0("~Condition", adjustment.formula)
            print(adjustment.formula)
            model <- data.frame(batchCombat, row.names = colnames(tabDF))
            design.mod.combat <- model.matrix(eval(parse(text = adjustment.formula)),
                                              data = model)
        }
        print(unique(batchCombat))
        batch_corr <- sva::ComBat(dat = tabDF, batch = batchCombat,
                                  mod = design.mod.combat, par.prior = TRUE, prior.plots = plot)
    }
    
    return(batch_corr)
}

# Function to run principal component analysis and plot correlations
runPCAandPlotCorrelationsRevised <- function(genesBySamples, samplesByCovariates, dataName, isKeyPlot=FALSE, 
                                             SCALE_DATA_FOR_PCA = TRUE, MIN_PVE_PCT_PC = 1.0, CORRELATION_TYPE = "pearson",
                                             ALSO_PLOT_ALL_COVARS_VS_PCA = TRUE, MAX_NUM_LEVELS_PER_COVAR = 50) {
    
    title = paste(ifelse(SCALE_DATA_FOR_PCA, "S", "Un-s"), "caled ", dataName, " ", " data in PCA; PVE >= ", MIN_PVE_PCT_PC, "%; ", CORRELATION_TYPE, " correlations ", sep="")
    writeLines(paste("\nRunning PCA and calculating correlations for:\n", title, sep=""))
    
    pcaRes <- runPCA(genesBySamples=genesBySamples, 
                     SCALE_DATA_FOR_PCA=SCALE_DATA_FOR_PCA, 
                     MIN_PVE_PCT_PC=MIN_PVE_PCT_PC)
    
    samplePCvals <- pcaRes$samplePCvals
    pve <- pcaRes$pve
    
    npca <- ncol(samplePCvals)
    
    colnames(samplePCvals) = paste(colnames(samplePCvals), " (", sprintf("%.2f", pve[1:npca]), "%)", sep="")
    
    # Find covariates without any missing data
    samplesByFullCovariates = samplesByCovariates[, which(apply(samplesByCovariates, 2, 
                                                                function(dat) all(!is.na(dat))))]
    EXCLUDE_VARS_FROM_FDR = setdiff(colnames(samplesByCovariates), colnames(samplesByFullCovariates))
    
    add_PC_res = list()
    significantCovars = c()
    
    LOOP_PLOT_ALL_COVARS = FALSE
    if (ALSO_PLOT_ALL_COVARS_VS_PCA) { LOOP_PLOT_ALL_COVARS = unique(c(LOOP_PLOT_ALL_COVARS, TRUE)) }
    
    for (PLOT_ALL_COVARS in LOOP_PLOT_ALL_COVARS) {
        corrRes = calcCompleteCorAndPlotRevised(samplePCvals, 
                                                samplesByCovariates, 
                                                CORRELATION_TYPE, 
                                                title, 
                                                WEIGHTS = pve[1:dim(samplePCvals)[2]],
                                                PLOT_ALL_COVARS, 
                                                EXCLUDE_VARS_FROM_FDR)
        add_PC_res[[length(add_PC_res)+1]] = list(plotData=corrRes$plot, isKeyPlot=(isKeyPlot && !PLOT_ALL_COVARS))
        if (!PLOT_ALL_COVARS) {
            significantCovars = corrRes$significantCovars
            Effects.significantCovars = corrRes$Effects.significantCovars
        }
    }
    
    return(list(significantCovars=significantCovars, PC_res=add_PC_res, Effects.significantCovars = Effects.significantCovars))
}

calcCompleteCorAndPlotRevised <- function(COMPARE_data,
                                          COVAR_data,
                                          correlationType, 
                                          title, 
                                          WEIGHTS = NULL,
                                          PLOT_ALL_COVARS=FALSE,
                                          EXCLUDE_VARS_FROM_FDR=NULL,
                                          MAX_FDR = 0.1) {
    
    # require(plyr)
    
    # Get factor and continuous covariates
    FactorCovariates <- colnames(COVAR_data)[sapply(COVAR_data,is.factor)]
    ContCovariates <- setdiff(colnames(COVAR_data),FactorCovariates)
    
    # Convert factor covariates to numeric vector
    # COVAR_data[,FactorCovariates] <- apply(COVAR_data[,FactorCovariates],2,
    #                                        function(x){x <- unclass(x)})
    
    for (var in FactorCovariates){
        COVAR_data[, var] <- unclass(COVAR_data[, var])
    }
    
    # Calculate correlation between compare_data and factor covariates
    if (length(FactorCovariates) > 0){
        comb <- expand.grid(colnames(COMPARE_data),FactorCovariates)
        factCont_cor <- apply(comb,1,
                              getFactorContAssociationStatistics,
                              cbind(COMPARE_data,COVAR_data[rownames(COMPARE_data),FactorCovariates]),
                              alpha=MAX_FDR)
        factCont_cor_vals <- matrix(factCont_cor['Estimate',],
                                    nrow = length(colnames(COMPARE_data)),
                                    ncol = length(FactorCovariates))
        factCont_cor_p <- matrix(factCont_cor['Pval',],
                                 nrow = length(colnames(COMPARE_data)),
                                 ncol = length(FactorCovariates))
        
        rownames(factCont_cor_vals) <- colnames(COMPARE_data)
        colnames(factCont_cor_vals) <- FactorCovariates
        
        rownames(factCont_cor_p) <- colnames(COMPARE_data)
        colnames(factCont_cor_p) <- FactorCovariates
    } else {
        factCont_cor_vals <- NULL
        factCont_cor_p <- NULL
    }
    
    # Calculate correlation between compare_data and factor covariates
    if (length(ContCovariates) > 0){
        cont_cor <- corr.test(COMPARE_data,
                              COVAR_data[,ContCovariates],
                              use='pairwise.complete.obs',
                              method=correlationType, 
                              adjust="none")
        cont_cor_vals <- cont_cor$r
        cont_cor_p <- cont_cor$p
        
        rownames(cont_cor_vals) <- colnames(COMPARE_data)
        colnames(cont_cor_vals) <- ContCovariates
        
        rownames(cont_cor_p) <- colnames(COMPARE_data)
        colnames(cont_cor_p) <- ContCovariates
    } else {
        cont_cor_vals <- NULL
        cont_cor_p <- NULL
    }
    
    all_cor_vals = cbind(factCont_cor_vals,cont_cor_vals)
    all_cor_p = cbind(factCont_cor_p,cont_cor_p)
    
    Effects.significantCovars = all_cor_vals
    Effects.significantCovars[all_cor_p>MAX_FDR] = 0
    Effects.significantCovars = colSums(abs(Effects.significantCovars)*replicate(dim(Effects.significantCovars)[2],WEIGHTS/sum(WEIGHTS)))
    Effects.significantCovars = Effects.significantCovars[order(abs(Effects.significantCovars),decreasing=T)]
    
    cor_mat = melt(all_cor_p, varnames=c("COMPARE", "COVAR"))
    colnames(cor_mat)[colnames(cor_mat) == "value"] = "pvalue"
    
    cor_mat$COMPARE = factor(cor_mat$COMPARE, levels=rownames(all_cor_p))
    cor_mat$COVAR = factor(cor_mat$COVAR, levels=colnames(all_cor_p))
    
    cor_mat$r = melt(all_cor_vals)$value
    
    calcFDRrows = rep(TRUE, nrow(cor_mat))
    markColumnsAsMissing = NULL
    if (!is.null(EXCLUDE_VARS_FROM_FDR)) {
        calcFDRrows = !(cor_mat$COVAR %in% EXCLUDE_VARS_FROM_FDR)
        markColumnsAsMissing = intersect(colnames(COVAR_data), EXCLUDE_VARS_FROM_FDR)
    }
    
    # Entries that pass the threshold of "significance":  
    markSignificantCorrelations = corMatFDRthreshFunc(cor_mat, indicesMask=calcFDRrows, MAX_FDR = 0.1)
    significantCorrelatedCovars = sort(unique(cor_mat$COVAR[markSignificantCorrelations]))
    
    markPotentialSignificantCorrelations = corMatFDRthreshFunc(cor_mat)
    # Specially mark only those incomplete covariates that would be significant in the context of all covariates:
    markPotentialSignificantCorrelations = markPotentialSignificantCorrelations & !calcFDRrows
    
    plotRows = 1:nrow(cor_mat)
    if (!PLOT_ALL_COVARS) {
        # Plot all correlations for:
        # 1) Covariates with at least one significant correlation
        # 2) Excluded covariates
        plotRows = (cor_mat$COVAR %in% significantCorrelatedCovars) | !calcFDRrows
    }
    plotCor = na.omit(cor_mat[plotRows, ])
    
    for (markCor in c("markSignificantCorrelations", "markPotentialSignificantCorrelations")) {
        useMarkCor = get(markCor)[plotRows]
        if (length(which(useMarkCor)) > 0) {
            plotCor[, markCor] = useMarkCor[ setdiff(1:length(useMarkCor), as.numeric(attr(plotCor, "na.action"))) ]
        }
    }
    
    if (!plyr::empty(plotCor)){
        plot = plotCorWithCompare(plotCor, title, paste("FDR <= ", MAX_FDR, sep=""), markColumnsAsMissing)
    } else{
        plot = NULL
    }
    
    return(list(plot=plot, significantCovars=as.character(significantCorrelatedCovars), Effects.significantCovars = Effects.significantCovars))
}













library(CovariateAnalysis)
library(reshape2)
library(ggplot2)
library(psych)
library(limma)
library(edgeR)
source('00-revised_functions.R')

# input.data <- readRDS('../Data/TCGA-LIHC-COUNTS.RDS')
# covariates <- read.csv('../Data/covariates.csv', row.names = 1)
# prime.variables <- c('plate', 'submitter_id', 'patient', 'shortLetterCode', 'definition', 'sample_type_id', 'sample_type', 'participant', 'condition')
# min.iter <- 20
# voom <- T
# out.dir <- '../Data'
# save.name <- 'test.RData'
# plot <- T

CovariatesCorrelationCorrection <- function(input.data, covariates, prime.variables, min.iter=20, voom=T, out.dir='', save.name=NULL, plot=F) {
    
    if (substring(out.dir, nchar(out.dir))!='/')
        out.dir <- paste0(out.dir, '/')
    if (!dir.exists(out.dir))
        stop('Output directory does not exist.')
    
    remove.variables <- c()
    
    for (col in colnames(covariates)){
        if (length(unique(na.omit(covariates[,col]))) == length(na.omit(covariates[,col])) || (length(unique(covariates[,col])) == 1)) {
            remove.variables <- c(remove.variables, col)      
        }
    }
    
    temp.covariates <- dplyr::select(covariates, -all_of(remove.variables))
    
    dataName <- ifelse(voom==T, 'NULL_design(voom-normalized)', 'NULL design')
    
    preAdjustedSigCovars <- runPCAandPlotCorrelationsRevised(input.data,
                                                             temp.covariates,
                                                             dataName, 
                                                             isKeyPlot=TRUE, 
                                                             MIN_PVE_PCT_PC = 1)
    
    if (plot==T){
        p <- preAdjustedSigCovars[["PC_res"]][[2]]$plotData
        ggsave(paste0(out.dir,"pre_correction_plot.pdf"),p, width = 18,height = 6)
    }
    
    residualCovars <- setdiff(preAdjustedSigCovars$significantCovars,prime.variables)
    residualSigCovars <- preAdjustedSigCovars
    covariatesEffects <- preAdjustedSigCovars$Effects.significantCovars[residualCovars]
    prime.variables <- unique(c(prime.variables, names(which.max(covariatesEffects))))
    
    count <- 0
    while(length(residualSigCovars$significantCovars)!=0 && count <= min.iter){
        writeLines(paste('Using following covariates in the model:',
                         paste(prime.variables, collapse=', '),
                         'as fixed effects'))
        
        DM <- getDesignMatrix(covariates[,prime.variables,drop=F],Intercept = T)
        DM$design = DM$design[,linColumnFinder(DM$design)$indepCols]
        
        if (voom==T){
            voom.weights <- voom(input.data, design=DM$design, plot=F)
        }
        else{
            voom.weights <- NULL
        }
        
        adjusted.fit <- lmFit(input.data,
                              design=DM$design,
                              weights=voom.weights)
        
        
        residual.exprs <- residuals.MArrayLM(adjusted.fit, input.data)
        
        residCovars <- setdiff(colnames(temp.covariates), prime.variables)
        
        temp.residual.exprs <- residual.exprs
        temp.residual.exprs[is.na(temp.residual.exprs)] = 0
        temp.residual.exprs <- temp.residual.exprs[apply(temp.residual.exprs,1, sd)!=0,]
        
        dataName <- ifelse(voom==T, 'adjusted design(voom-normalized)', 'adjusted design')
        
        residualSigCovars = runPCAandPlotCorrelationsRevised(temp.residual.exprs, 
                                                             temp.covariates[, residCovars, drop=F], 
                                                             dataName,
                                                             isKeyPlot=TRUE)
        
        residCovars = setdiff(residualSigCovars$significantCovars, prime.variables)
        covariatesEffects = residualSigCovars$Effects.significantCovars[residCovars]
        prime.variables <- unique(c(prime.variables, names(which.max(covariatesEffects))))
    }
    
    if (plot==T){
        residualSigCovars = runPCAandPlotCorrelationsRevised(temp.residual.exprs, 
                                                             temp.covariates,
                                                             dataName,
                                                             isKeyPlot=TRUE)
        
        p <- residualSigCovars[["PC_res"]][[2]]$plotData
        
        ggsave(paste0(out.dir,"post_correction_plot.pdf"),p,width = 18,height = 6)
    }
    
    DM <- getDesignMatrix(covariates[, prime.variables,drop = F], Intercept = F)
    
    DM$design <- DM$design[,linColumnFinder(DM$design)$indepCols]
    
    FIT <- lmFit(input.data, DM$design)
    
    residual.exprs <- residuals.MArrayLM(FIT,input.data)
    
    output <- list('primeVariables'=prime.variables, 'residualExprs'=residual.exprs)
    
    if(!is.null(save.name)){
        save(output, paste0(out.dir, save.name))
    }
    
    return (output)
}

