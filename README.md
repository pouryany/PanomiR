<!-- badges: start -->

[![R-CMD-check](https://github.com/pouryany/PanomiR/workflows/R-CMD-check/badge.svg)](https://github.com/pouryany/PanomiR/actions)
[![lint](https://github.com/pouryany/PanomiR/workflows/lint/badge.svg)](https://github.com/pouryany/PanomiR/actions)
<!-- badges: end -->

## Introduction

PanomiR is a package for pathway and microRNA Analysis of gene
expression data. This document provides details about how to install and
utilize various functionality in PanomiR.

For questions, comments, and other queries, contact <pouryany@gmail.com>

## Installation

PanomiR can be accessed via Bioconductor. To install, start R and run
the following code.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PanomiR")
```

You can also install the latest development version of PanomiR using
GitHub.

``` r
devtools::install_github("pouryany/PanomiR")
```

## Overview

PanomiR is a pipeline to prioritize disease-associated miRNAs based on
activity of disease-associated pathways. The input datasets for PanomiR
are (a) a gene expression disease dataset along with covariates, (b) a
background collection of pathways/genesets, and (c) a collection of
miRNAs containing gene targets.

The general workflow of PanomiR is (a) generation of pathway summary
statistics from gene expression data, (b) detection of differentially
activated pathways, (c) finding coherent groups, or clusters, of
differentially activated pathways, and (d) detecting miRNAs targeting
each group of pathways.

Individual steps of the workflow can be used in isolation to carry out
different analyses. The following sections outline each step and
material needed to execute PanomiR.

## 1. Pathway summarization

PanomiR can generate pathway activity profiles given a gene expression
dataset and a list of pathways.Pathway summaries are numbers that
represent the overall activity of genes that belong to each pathway.
These numbers are calculated based on a methodology previously described
in part in (Altschuler et al. 2013; Joachim et al. 2018). Briefly, genes
in each sample are ranked by their expression values and then pathway
summaries are calculated as the average rank-squared of genes within a
pathway. The summaries are then center and scaled (zNormalized) across
samples.

The default list of background pathways in PanomiR is formatted into a
table (`data("path_gene_table")`). The table is based on canonical
pathways collection of Molecular Signatures Database (MSigDB) V6.2 and
it contains annotated pathways from a variety of sources (Liberzon et
al. 2011).

\*\* Users interested in using other pathway/geneset backgrounds, such
as newer versions of MSigDB or KEGG, should refer to the
[appendix](#geneset) of this manual.

This section uses a reduced example dataset from The Cancer Genome Atlas
(TCGA) Liver Hepatocellular Carcinoma (LIHC) dataset to generate Pathway
summary statistics (Ally et al. 2017). **Note:** Make sure that you
select gene representation type that matches the rownames of your
expression data. The type can be modified using the `id` argument in the
function below. The default value for this argument is `ENSEMBL`.

``` r
library(PanomiR)

# Pathway reference from the PanomiR package
data("path_gene_table")
data("miniTestsPanomiR")
# Generating pathway summary statistics 

summaries <- pathwaySummary(miniTestsPanomiR$mini_LIHC_Exp,
                            path_gene_table, method = "x2",
                            zNormalize = TRUE, id = "ENSEMBL")

head(summaries)[,1:2]
```

    ##                         TCGA-BC-A10S-01A-22R-A131-07
    ## BIOCARTA_41BB_PATHWAY                     -0.1506216
    ## BIOCARTA_ACE2_PATHWAY                     -0.5676447
    ## BIOCARTA_ACH_PATHWAY                      -0.3211747
    ## BIOCARTA_ACTINY_PATHWAY                    1.4363526
    ## BIOCARTA_AGPCR_PATHWAY                    -0.1948523
    ## BIOCARTA_AGR_PATHWAY                       0.6802993
    ##                         TCGA-BC-4073-01B-02R-A131-07
    ## BIOCARTA_41BB_PATHWAY                     -0.1269436
    ## BIOCARTA_ACE2_PATHWAY                     -0.8327436
    ## BIOCARTA_ACH_PATHWAY                      -0.4390042
    ## BIOCARTA_ACTINY_PATHWAY                    1.4975456
    ## BIOCARTA_AGPCR_PATHWAY                    -0.2499193
    ## BIOCARTA_AGR_PATHWAY                       0.5420588

## 2. Differential Pathway activation

Once you generate the pathway activity profiles, as discussed in the
last section, there are several analysis that you can perform. We have
bundled some of the most important ones into standalone functions. Here,
we describe differential pathway activation profiling, which is
examining differences in pathway activity profiles in user-determined
conditions.

At this stage you need to provide a pathway-gene association table, an
expression dataset, and a covariates table. You need to specity what
covariates you would like to contrast. You also need to provide a
contrast, as formatted in limma. If the contrast is not provided, the
function assumes the first two levels of the provided contrast
covariate. **Note:** make sure the contrast covariate is formatted as
factor.

``` r
output0 <- differentialPathwayAnalysis(
                        geneCounts = miniTestsPanomiR$mini_LIHC_Exp,
                        pathways =  path_gene_table,
                        covariates = miniTestsPanomiR$mini_LIHC_Cov,
                        condition = 'shortLetterCode')

de.paths <- output0$DEP

head(de.paths,3)
```

    ##                                                 logFC   AveExpr          t
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING -0.9159376 0.3044281 -10.404966
    ## BIOCARTA_AKT_PATHWAY                       -0.5744103 0.3123897  -6.770069
    ## PID_IL5_PATHWAY                            -0.6219876 0.4240432  -6.255756
    ##                                                 P.Value   adj.P.Val        B
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING 1.942463e-06 0.002012391 5.240095
    ## BIOCARTA_AKT_PATHWAY                       6.903010e-05 0.035757593 2.126311
    ## PID_IL5_PATHWAY                            1.276971e-04 0.040289104 1.550780
    ##                                                                       contrast
    ## REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING shortLetterCodeTP-shortLetterCodeNT
    ## BIOCARTA_AKT_PATHWAY                       shortLetterCodeTP-shortLetterCodeNT
    ## PID_IL5_PATHWAY                            shortLetterCodeTP-shortLetterCodeNT

## 3. Finding clusters of pathways

PanomiR provides a function to find groups coordinated differentially
activated pathways based on a pathway co-expression network (PCxN)
previously described in (Pita-Juárez et al. 2018). Briefly, PCxN is a
network where nodes are pathways and links are co-expression between the
nodes. It is formatted into a table were rows represent edges. The edges
of PCxN are marked by two numbers, 1- a correlation co-efficient and 2-
a significance adjusted p-value. Cut-offs for both of these numbers can
be manually set using PanomiR functions. See function manuals for more
info.

PCxN and its associated genesets are already released and can be
accessed through following Bioconductor packages:
[pcxn](http://bioconductor.org/packages/release/bioc/html/pcxn.html) and
[pcxnData](http://bioconductor.org/packages/release/data/experiment/html/pcxnData.html).

Here we have provided a small version of PCxN for tutorial purposes. A
more recent version of PCxN based on MSigDB V6.2 is available through
the data repository accompanying PanomiR manuscript, which can be found
[here](https://github.com/pouryany/PanomiR_paper).

``` r
# using an updated version of pcxn 

set.seed(2)
pathwayClustsLIHC <- mappingPathwaysClusters(
                            pcxn = miniTestsPanomiR$miniPCXN, 
                            dePathways = de.paths[1:300,],
                            topPathways = 200,
                            outDir=".",
                            plot = FALSE,
                            subplot = FALSE,
                            prefix='',
                            clusteringFunction = "cluster_louvain",
                            correlationCutOff = 0.1)


head(pathwayClustsLIHC$Clustering)
```

    ##                      Pathway cluster
    ## 1       BIOCARTA_NO1_PATHWAY       1
    ## 2       BIOCARTA_AKT_PATHWAY       1
    ## 3       BIOCARTA_ALK_PATHWAY       1
    ## 4     BIOCARTA_RANKL_PATHWAY       1
    ## 5       BIOCARTA_MCM_PATHWAY       3
    ## 6 BIOCARTA_CELLCYCLE_PATHWAY       3

## 4. Prioritizing miRNAs per cluster of pathways.

PanomiR identifies miRNAs that target clusters of pathways, as defined
in the last section. In order to this, you would need a reference table
of miRNA-Pathway association score (enrichment). We recommend using a
customized miRNA-Pathway association table, tailored to your
experimental data. This section provides an overview of prioritization
process. Readers interested in knowing more about the technical details
of PanomiR are refered to accompaniying publication (Work under
preparation).

### Enrichment reference

Here, we provide a preprocessed small example table of miRNA-pathway
enrichment in `miniTestsPanomiR$miniEnrich` object. This table contains
enrichment analysis results using Fisher’s Exact Test between MSigDB
pathways and TargetScan miRNA targets. The individual components are
accessible via `data(msigdb_c2)` and `data(targetScan_03)` (Agarwal et
al. 2015; Liberzon et al. 2011). This example table is contains only a
full subset of the full pairwise enrichment. You can refer to [section
5](#geneset) of this manual on how to create full tables and how to
customize them to your specific gene expression data.

### Generating targeting scores

PanomiR generates a score for individual miRNAs targeting a group of
pathways. These scores are generated based on the reference enrichment
table. We are interested in knowing to what extent each miRNA targets
pathway clusters identified in the last step (see previous section).
PanomiR constructs a null distribution of this targeting score for each
miRNA. The significance of observed scores from a given group of
pathways (clusters in this case) is contrasted against the null
distribution to generate a targeting p-value. These p-values are used to
rank miRNAs per cluster.

### Sampling parameter

The above described process requires repeated sampling to empirically
obtain the null distribution. The argument `sampRate` denotes the number
of repeats in the process. Note that in the example below, we use a
sampling rate of 50, the recommended rate is between 500-1000. Also, we
set the saveSampling argument to FALSE. This argument, if set TRUE,
ensures that the null distribution is obtain only once. This argument
should be set to TRUE if you wish to save your sampling and check for
different outputs from the clustering algorithms or pathway thresholds.

``` r
set.seed(1)
output2 <- prioritizeMicroRNA(enriches0 = miniTestsPanomiR$miniEnrich,
                    pathClust = miniTestsPanomiR$miniPathClusts$Clustering,
                    topClust = 1,
                    sampRate = 50, 
                    method = c("aggInv"),
                    outDir = "Output/",
                    dataDir = "outData/",
                    saveSampling = FALSE,
                    runJackKnife = FALSE,
                    numCores = 1,
                    prefix = "outmiR",
                    saveCSV = FALSE)
```

    ## Working on Cluster1.

    ## Performing aggInv function.

    ## aggInv Method Done

``` r
head(output2$Cluster1)
```

    ##                                 x cluster_hits aggInv_cover  aggInv_pval
    ## 1                hsa-miR-101-3p.2            6   -1.9566603 0.0001216703
    ## 2                hsa-miR-101-3p.1            4   -0.3395771 0.0006214715
    ## 3 hsa-miR-124-3p.2/hsa-miR-506-3p            7   -0.2357761 0.0008599272
    ## 4                 hsa-miR-1247-5p            4   -1.6599230 0.0021625662
    ## 5                 hsa-miR-1249-3p            1   -2.4578993 0.0042061415
    ## 6                 hsa-miR-1252-5p            4   -0.7572036 0.0050836835
    ##    aggInv_fdr
    ## 1 0.002433406
    ## 2 0.005732848
    ## 3 0.005732848
    ## 4 0.010812831
    ## 5 0.016824566
    ## 6 0.016945612

## 5. miRNA-Pathway enrichment tables

PanomiR best performs on tissue/experiment-customized datasets. In order
to do this, you need to create a customized enrichment table. You can
simply do so by using the pathway and miRNA list that we have provided
as a part of the package. simply, plug in the name of the genes present
(expressed) in your experiment in the following code

``` r
# using an updated version of pcxn 
data("msigdb_c2")
data("targetScan_03")


customeTableEnrich <- miRNAPathwayEnrichment(mirSets = targetScan_03,
                                              pathwaySets = msigdb_c2,
                                              geneSelection = yourGenes,
                                              mirSelection = yourMicroRNAs,
                                              fromID = "ENSEMBL",
                                              toID = "ENTREZID",
                                              minPathSize = 9,
                                              numCores = 1,
                                              outDir = ".",
                                              saveOutName = NULL)
```

In the above section, the field `fromID` denotes the gene representation
format of your input list. Here is a quick example that runs fast. Note
that the `miRNAPathwayEnrichment()` function creates a detailed report
with parameters that are used internally. To get a smaller table that is
suitable for publication purposes, use `reportEnrichment()` function.

``` r
# using an updated version of pcxn 
data("msigdb_c2")
data("targetScan_03")

tempEnrich <-miRNAPathwayEnrichment(targetScan_03[1:30],msigdb_c2[1:30])

head(reportEnrichment(tempEnrich))
```

    ##                               miRNA                Pathway    pval pAdjust
    ## 14                  hsa-miR-1252-5p   BIOCARTA_CSK_PATHWAY 0.00179   0.259
    ## 55    hsa-miR-1271-5p/hsa-miR-96-5p   BIOCARTA_AKT_PATHWAY 0.01100   1.000
    ## 122                hsa-miR-124-3p.1 BIOCARTA_RANKL_PATHWAY 0.01550   1.000
    ## 53  hsa-miR-124-3p.2/hsa-miR-506-3p   BIOCARTA_AKT_PATHWAY 0.03740   1.000
    ## 99                  hsa-miR-1252-5p BIOCARTA_AGPCR_PATHWAY 0.03940   1.000
    ## 112                hsa-miR-124-3p.1   BIOCARTA_BCR_PATHWAY 0.05360   1.000

## 6. Customized genesets and recommendations

PanomiR can integrate genesets and pathways from external sources
including those annotated in MSigDB. In order to do so, you need to
provide a `GeneSetCollection` object as defined in the `GSEABase`
package.

The example below illustrates how to use external sources to create your
own customized pathway-gene association table. This customized can then
replaced `path_gene_table` input in functions described in sections 1,2,
and 5 of this manual.

``` r
data("gscExample")

newPathGeneTable <-tableFromGSC(gscExample)
```

    ## 

    ## 

    ## 'select()' returned 1:1 mapping between keys and columns
    ## 'select()' returned 1:1 mapping between keys and columns

``` r
head(newPathGeneTable)
```

    ## # A tibble: 6 × 3
    ##   Pathway                               ENTREZID ENSEMBL        
    ##   <chr>                                 <chr>    <chr>          
    ## 1 NIKOLSKY_BREAST_CANCER_16Q24_AMPLICON 124245   ENSG00000158545
    ## 2 NIKOLSKY_BREAST_CANCER_16Q24_AMPLICON 1800     ENSG00000015413
    ## 3 NIKOLSKY_BREAST_CANCER_16Q24_AMPLICON 404550   ENSG00000154102
    ## 4 NIKOLSKY_BREAST_CANCER_16Q24_AMPLICON 6687     ENSG00000197912
    ## 5 NIKOLSKY_BREAST_CANCER_16Q24_AMPLICON 750      ENSG00000221819
    ## 6 NIKOLSKY_BREAST_CANCER_16Q24_AMPLICON 2588     ENSG00000141012

The the pathway correlation network in section 3 is build upon an MSigDB
V6.2, canonical pathways (cp) collection dataset that includes KEGG
Pathways. KEGG prohibits distribution of its pathways by third parties.
However, use can access desired versions of MSigDB in gmt format via
[this link](https://www.gsea-msigdb.org/gsea/downloads_archive.jsp)
(Subramanian et al. 2005).

The library `msigdb` provides an programmatic interface to download
different geneset collections. Including how to add KEGG pathways or
download mouse genesets. Use the this [MSigDB
tutorial](https://bioconductor.org/packages/release/data/experiment/vignettes/msigdb/inst/doc/msigdb.html)
to create your desired gene sets.

You can also use the following code chunk to create pathway-gene
association tables from gmt files.

``` r
library(GSEABase)

yourGeneSetCollection <- getGmt("YOUR GMT FILE")
newPathGeneTable      <- tableFromGSC(yourGeneSetCollection)
```

## Session info

``` r
sessionInfo()
```

    ## R Under development (unstable) (2021-12-03 r81290)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur/Monterey 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] PanomiR_0.99.4
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fgsea_1.21.0           colorspace_2.0-2       ggtree_3.3.1          
    ##   [4] ellipsis_0.3.2         qvalue_2.27.0          XVector_0.35.0        
    ##   [7] aplot_0.1.1            rstudioapi_0.13        farver_2.1.0          
    ##  [10] graphlayouts_0.8.0     ggrepel_0.9.1          bit64_4.0.5           
    ##  [13] scatterpie_0.1.7       AnnotationDbi_1.57.1   fansi_0.5.0           
    ##  [16] splines_4.2.0          cachem_1.0.6           GOSemSim_2.21.0       
    ##  [19] knitr_1.37             polyclip_1.10-0        jsonlite_1.7.2        
    ##  [22] annotate_1.73.0        GO.db_3.14.0           png_0.1-7             
    ##  [25] graph_1.73.0           ggforce_0.3.3          compiler_4.2.0        
    ##  [28] httr_1.4.2             assertthat_0.2.1       Matrix_1.4-0          
    ##  [31] fastmap_1.1.0          lazyeval_0.2.2         cli_3.1.0             
    ##  [34] limma_3.51.2           tweenr_1.0.2           htmltools_0.5.2       
    ##  [37] tools_4.2.0            igraph_1.2.11          gtable_0.3.0          
    ##  [40] glue_1.6.0             GenomeInfoDbData_1.2.7 reshape2_1.4.4        
    ##  [43] DO.db_2.9              dplyr_1.0.7            fastmatch_1.1-3       
    ##  [46] Rcpp_1.0.7             enrichplot_1.15.2      Biobase_2.55.0        
    ##  [49] vctrs_0.3.8            Biostrings_2.63.1      ape_5.6-1             
    ##  [52] nlme_3.1-153           ggraph_2.0.5           xfun_0.29             
    ##  [55] stringr_1.4.0          lifecycle_1.0.1        clusterProfiler_4.3.1 
    ##  [58] XML_3.99-0.8           DOSE_3.21.2            org.Hs.eg.db_3.14.0   
    ##  [61] zlibbioc_1.41.0        MASS_7.3-54            scales_1.1.1          
    ##  [64] tidygraph_1.2.0        parallel_4.2.0         RColorBrewer_1.1-2    
    ##  [67] yaml_2.2.1             memoise_2.0.1          gridExtra_2.3         
    ##  [70] ggplot2_3.3.5          downloader_0.4         ggfun_0.0.4           
    ##  [73] yulab.utils_0.0.4      stringi_1.7.6          RSQLite_2.2.9         
    ##  [76] S4Vectors_0.33.9       tidytree_0.3.6         BiocGenerics_0.41.2   
    ##  [79] BiocParallel_1.29.10   GenomeInfoDb_1.31.1    rlang_0.4.12          
    ##  [82] pkgconfig_2.0.3        bitops_1.0-7           evaluate_0.14         
    ##  [85] lattice_0.20-45        purrr_0.3.4            treeio_1.19.1         
    ##  [88] patchwork_1.1.1        shadowtext_0.1.0       bit_4.0.4             
    ##  [91] tidyselect_1.1.1       GSEABase_1.57.0        plyr_1.8.6            
    ##  [94] magrittr_2.0.1         R6_2.5.1               IRanges_2.29.1        
    ##  [97] generics_0.1.1         DBI_1.1.2              pillar_1.6.4          
    ## [100] withr_2.4.3            KEGGREST_1.35.0        RCurl_1.98-1.5        
    ## [103] tibble_3.1.6           crayon_1.4.2           utf8_1.2.2            
    ## [106] rmarkdown_2.11         viridis_0.6.2          grid_4.2.0            
    ## [109] data.table_1.14.2      blob_1.2.2             forcats_0.5.1         
    ## [112] digest_0.6.29          xtable_1.8-4           tidyr_1.1.4           
    ## [115] gridGraphics_0.5-1     stats4_4.2.0           munsell_0.5.0         
    ## [118] viridisLite_0.4.0      ggplotify_0.1.0

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-agarwal2015predicting" class="csl-entry">

Agarwal, Vikram, George W Bell, Jin-Wu Nam, and David P Bartel. 2015.
“Predicting Effective microRNA Target Sites in Mammalian mRNAs.” *Elife*
4: e05005.

</div>

<div id="ref-ally2017comprehensive" class="csl-entry">

Ally, Adrian, Miruna Balasundaram, Rebecca Carlsen, Eric Chuah, Amanda
Clarke, Noreen Dhalla, Robert A Holt, et al. 2017. “Comprehensive and
Integrative Genomic Characterization of Hepatocellular Carcinoma.”
*Cell* 169 (7): 1327–41.

</div>

<div id="ref-altschuler2013pathprinting" class="csl-entry">

Altschuler, Gabriel M, Oliver Hofmann, Irina Kalatskaya, Rebecca Payne,
Shannan J Ho Sui, Uma Saxena, Andrei V Krivtsov, et al. 2013.
“Pathprinting: An Integrative Approach to Understand the Functional
Basis of Disease.” *Genome Medicine* 5 (7): 1–13.

</div>

<div id="ref-joachim2018relative" class="csl-entry">

Joachim, Rose B, Gabriel M Altschuler, John N Hutchinson, Hector R Wong,
Winston A Hide, and Lester Kobzik. 2018. “The Relative Resistance of
Children to Sepsis Mortality: From Pathways to Drug Candidates.”
*Molecular Systems Biology* 14 (5): e7998.

</div>

<div id="ref-liberzon2011molecular" class="csl-entry">

Liberzon, Arthur, Aravind Subramanian, Reid Pinchback, Helga
Thorvaldsdóttir, Pablo Tamayo, and Jill P Mesirov. 2011. “Molecular
Signatures Database (MSigDB) 3.0.” *Bioinformatics* 27 (12): 1739–40.

</div>

<div id="ref-pita2018pathway" class="csl-entry">

Pita-Juárez, Yered, Gabriel Altschuler, Sokratis Kariotis, Wenbin Wei,
Katjuša Koler, Claire Green, Rudolph E Tanzi, and Winston Hide. 2018.
“The Pathway Coexpression Network: Revealing Pathway Relationships.”
*PLoS Computational Biology* 14 (3): e1006042.

</div>

<div id="ref-subramanian2005gene" class="csl-entry">

Subramanian, Aravind, Pablo Tamayo, Vamsi K Mootha, Sayan Mukherjee,
Benjamin L Ebert, Michael A Gillette, Amanda Paulovich, et al. 2005.
“Gene Set Enrichment Analysis: A Knowledge-Based Approach for
Interpreting Genome-Wide Expression Profiles.” *Proceedings of the
National Academy of Sciences* 102 (43): 15545–50.

</div>

</div>
