
library(dplyr)
enriches <- readRDS("../test_cases/LIHC/Data/LIHCGenesLIHCMirsENRICHMENT_Tarbase.RDS")



enriches <- unique(enriches[,c("x","mirset_Size")])

enriches0 <- read.csv("../test_cases/LIHC/Output/cluster_fast_greedy_Prioritization_Tarbase/x2_LIHCGene_Tarbase1000_samples_clustNo_3.csv")


enriches2  <- left_join(enriches0,enriches)


plot(log(enriches2$mirset_Size),enriches2$cluster_hits,)

cor1 <- cor.test(log2(enriches3$mirset_Size),-log2(enriches3$AggInv_pval),method = "spearman")
cor2 <- cor.test(log2(enriches2$mirset_Size),enriches2$cluster_hits,method = "spearman")


library(ggplot2)
library(cowplot)


enriches3 <- enriches2#[enriches2$AggInv_fdr <1,]
cor1      <- cor.test(log2(enriches3$mirset_Size),
                      -log2(enriches3$AggInv_pval),
                      method = "spearman")

ggplot(enriches3,aes(log2(mirset_Size),-log2(AggInv_pval))) +
    geom_point(alpha=0.8) + 
    theme_cowplot() +
    xlab("log(Number of miRNA target)")+
    ylab(paste0( "- log(caraway p-value)")) +
    labs(title = paste0("CARAWAY Pathway Size vs P-value"))+
    theme(plot.title = element_text(size=22),
          axis.text=element_text(size=14),
          axis.text.x = element_text(face = "bold",size = 14, angle = 0,
                                     vjust = 1,hjust = 1),
          axis.text.y =  element_text(size = 20),
          axis.title=element_text(size=14,face="bold"),
          strip.text = element_text(size = 14,angle = 90),
          legend.position = "none") +
     annotate("text", x = 10, y = 100, 
              label = paste0("italic(Correlation) == ", signif(cor1$estimate,2) ),
              parse = T,
              size  = 7)


enriches3 <- enriches2
cor2      <- cor.test(log2(enriches3$mirset_Size),
                      (enriches3$cluster_hits),
                      method = "spearman")



ggplot(enriches3,aes(log2(mirset_Size),(cluster_hits))) +
    geom_point(alpha=0.6) + 
    theme_cowplot() +
    xlab("log(Number of miRNA target)")+
    ylab(paste0( "Number of enriched Pathways")) +
    labs(title = paste0("Enrichment Pathway Size vs P-value"))+
    theme(plot.title = element_text(size=22),
          axis.text=element_text(size=14),
          axis.text.x = element_text(face = "bold",size = 14, angle = 0,
                                     vjust = 1,hjust = 1),
          axis.text.y =  element_text(size = 20),
          axis.title=element_text(size=14,face="bold"),
          strip.text = element_text(size = 14,angle = 90),
          legend.position = "none")+
    annotate("text", x = 5, y = 40, 
             label = paste0("italic(Correlation) == ", signif(cor2$estimate,2) ),
             parse = T,
             size  = 7)




