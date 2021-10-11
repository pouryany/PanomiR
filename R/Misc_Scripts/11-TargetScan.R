rm(list = ls())
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)


targetScan0 <- data.table::fread("../../Data/Summary_Counts.all_predictions.txt")
filterScan  <- grep("9606",targetScan0$`Species ID`)
targetScan  <- targetScan0[filterScan,]
#targetScan0 <- targetScan0[filterScan,]



family_info <- data.table::fread("../../Data/TargetScan72_July2021/miR_Family_Info.txt") 
family_info <- family_info[grep("9606",family_info$`Species ID`),]
family_info <- family_info[,-c("Species ID")]
targetScan  <- left_join(targetScan,family_info,by = c("Representative miRNA" = "MiRBase ID"))


targetScan$`Cumulative weighted context++ score` <- as.numeric(targetScan$`Cumulative weighted context++ score`)

targetScan  <- targetScan[!is.na(targetScan$`Cumulative weighted context++ score`),]
targetScan  <- targetScan[(targetScan$`Cumulative weighted context++ score`) <0,]

targetScan  <- targetScan[targetScan$`Family Conservation?` != "-1",]


target.con  <- targetScan[targetScan$`Family Conservation?` != "0",]
target.con  <- target.con[target.con$`Total num conserved sites` !=0,]
target.poor <- targetScan[targetScan$`Family Conservation?` == "0",]

target.poor <- target.poor[target.poor$`Cumulative weighted context++ score` <= -0.20,]
#targetScan <- targetScan[targetScan$`Total num conserved sites` !=0,]
targetScan0 <- rbind(target.con,target.poor)

targetScan1 <- targetScan[(targetScan$`Cumulative weighted context++ score`) < -0.1,]
targetScan2 <- targetScan[(targetScan$`Cumulative weighted context++ score`) < -0.2,]
targetScan3 <- targetScan[(targetScan$`Cumulative weighted context++ score`) < -0.3,]

scenarios <- list("00" = targetScan0,
                  "01" = targetScan1,
                  "02" = targetScan2,
                  "03" = targetScan3)


targetSelector  <- function(targetScan1){
    ttt <- targetScan1 %>% group_by(.,`miRNA family`) %>%
        dplyr::select(.,`Representative miRNA`) %>%
        unique() %>%
        mutate(superFamily = paste0(`Representative miRNA`,collapse = "/")) 
    
    ttt <- ttt[!duplicated(ttt$`miRNA family`),]
    ttt <- ttt[,c(1,3)]
    
    ttt <- dplyr::left_join(targetScan1,ttt)
    
    ttt       <- ttt %>% dplyr::select(.,`Gene Symbol`,superFamily) %>%
                    group_split(.,superFamily) 
    mir.list  <- sapply(ttt, function(X){dplyr::slice(X,1) %>% 
                    dplyr::select(.,superFamily) %>% pull})
    ttt       <- sapply(ttt,function(X){dplyr::select(X,`Gene Symbol`) %>%
                                            unique})
    
    names(ttt) <- mir.list
    
    mir.sets2  <- lapply(ttt, function(X){
        gene.df <- bitr(X, fromType ="SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)
        return(gene.df$ENTREZID)
    })
    
    return(mir.sets2)
}



for(i in 1:length(scenarios)){
    mir.sets2  <- targetSelector(scenarios[[i]])
    
    saveRDS(mir.sets2,paste0("../../Data/preprocessed",
                             "/NORMALIZED_MIRSETS_TargetScan",
                             names(scenarios)[i],
                             ".rds"))
}

