# libraries needed
library(dplyr)
# library(plyr)
library(caret)
library(psych)
library(pROC)

SPLIT = "4"
thres <- 0.1
ressave <- T
ressuf <- "_v6"

TRAIN_SPLIT = paste0("train",SPLIT)
TEST_SPLIT = paste0("test",SPLIT)

#load data
#################################### Load Data ##########################################################

load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/MainModel_phenotypes_893people_12pheno_",TRAIN_SPLIT,"_v1.RData"))
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/ExtendedModel_phenotypes_1117people_48pheno_v2.RData")
extendedmodel_train <- extendedmodel[rownames(mainmodel_train),]
load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_strains_893people_",TRAIN_SPLIT,"_v1.RData"))
load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_893people_",TRAIN_SPLIT,"_v1.RData"))
load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKO_893people_",TRAIN_SPLIT,"_v1.RData"))
load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKEGGModule_893people_",TRAIN_SPLIT,"_v1.RData"))

load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/MainModel_phenotypes_224people_12pheno_",TEST_SPLIT,"_v1.RData"))
extendedmodel_test <- extendedmodel[rownames(mainmodel_test),]
load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_strains_224people_",TEST_SPLIT,"_v1.RData"))
load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_224people_",TEST_SPLIT,"_v1.RData"))
load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKO_224people_",TEST_SPLIT,"_v1.RData"))
load(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKEGGModule_224people_",TEST_SPLIT,"_v1.RData"))

mainmodel <- mainmodel_train
extendedmodel <- extendedmodel_train
microbes <- microbes_train
pathways <- pathways_train
pathwaysKO <- pathwaysKO_train
pathwaysKEGGModule <- pathwaysKEGGModule_train


################################### train ###########################################################
#merge bacteria (species) to genus, family
cn0 <- colnames(microbes)
cn <- strsplit(cn0,";",fixed = T)
cn <- as.data.frame(cn); colnames(cn) <- cn0; cn <- as.data.frame(t(cn))
colnames(cn) <- c("kingdom","phylum","class","order","family","genus","species")

microbes_toagg0 <- as.data.frame(t(microbes))
microbes_toagg0 <- cbind(cn,microbes_toagg0);
# microbes_toagg0 <- microbes_toagg0[microbes_toagg0$genus=="g__Lawsonia",]

#genus
microbes_toagg <- as.data.frame(t(microbes))
stopifnot(all(rownames(microbes_toagg)==rownames(cn)))
microbes_toagg <- cbind(cn$genus,microbes_toagg); cnn <- colnames(microbes_toagg); cnn[1] <- "toAgg"; colnames(microbes_toagg) <- cnn
microbes_toagg <- microbes_toagg[microbes_toagg$toAgg!="g__",]

microbes_genus <- microbes_toagg %>% 
  group_by(toAgg) %>% 
  summarise(across(everything(), ~sum(.x,na.rm = TRUE)))
microbes_genus <- as.data.frame(microbes_genus)
rownames(microbes_genus) <- microbes_genus$toAgg
microbes_genus <- microbes_genus[,2:length(microbes_genus)]
microbes_genus <- as.data.frame(t(microbes_genus))
expcn <- colnames(microbes_genus); expcn<-gsub("g__","Genus__",expcn,fixed = T);colnames(microbes_genus)<-expcn

#family
microbes_toagg <- as.data.frame(t(microbes))
stopifnot(all(rownames(microbes_toagg)==rownames(cn)))
microbes_toagg <- cbind(cn$family,microbes_toagg); cnn <- colnames(microbes_toagg); cnn[1] <- "toAgg"; colnames(microbes_toagg) <- cnn
microbes_toagg <- microbes_toagg[microbes_toagg$toAgg!="f__",]

microbes_family <- microbes_toagg %>% 
  group_by(toAgg) %>% 
  summarise(across(everything(), ~sum(.x,na.rm = TRUE)))
microbes_family <- as.data.frame(microbes_family)
rownames(microbes_family) <- microbes_family$toAgg
microbes_family <- microbes_family[,2:length(microbes_family)]
microbes_family <- as.data.frame(t(microbes_family))
expcn <- colnames(microbes_family); expcn<-gsub("f__","Family__",expcn,fixed = T);colnames(microbes_family)<-expcn

##############################################################################################
# convert to relative abundance
microbes_species <- t(microbes)
microbes_species <- scale(microbes_species,center = F,scale = colSums(microbes_species,na.rm = T))
microbes_species <- as.data.frame(t(microbes_species))

microbes_genus <- t(microbes_genus)
microbes_genus <- scale(microbes_genus,center = F,scale = colSums(microbes_genus,na.rm = T))
microbes_genus <- as.data.frame(t(microbes_genus))

microbes_family <- t(microbes_family)
microbes_family <- scale(microbes_family,center = F,scale = colSums(microbes_family))
microbes_family <- as.data.frame(t(microbes_family))

microbes_combined <- cbind(microbes_species,microbes_genus,microbes_family)
stopifnot(all(rownames(microbes_combined)==rownames(extendedmodel)))

################################### test ###########################################################
microbes <- microbes_test

#merge bacteria (species) to genus, family
cn0 <- colnames(microbes)
cn <- strsplit(cn0,";",fixed = T)
cn <- as.data.frame(cn); colnames(cn) <- cn0; cn <- as.data.frame(t(cn))
colnames(cn) <- c("kingdom","phylum","class","order","family","genus","species")

microbes_toagg0 <- as.data.frame(t(microbes))
microbes_toagg0 <- cbind(cn,microbes_toagg0);
# microbes_toagg0 <- microbes_toagg0[microbes_toagg0$genus=="g__Lawsonia",]

#genus
microbes_toagg <- as.data.frame(t(microbes))
stopifnot(all(rownames(microbes_toagg)==rownames(cn)))
microbes_toagg <- cbind(cn$genus,microbes_toagg); cnn <- colnames(microbes_toagg); cnn[1] <- "toAgg"; colnames(microbes_toagg) <- cnn
microbes_toagg <- microbes_toagg[microbes_toagg$toAgg!="g__",]

microbes_genus <- microbes_toagg %>% 
  group_by(toAgg) %>% 
  summarise(across(everything(), ~sum(.x,na.rm = TRUE)))
microbes_genus <- as.data.frame(microbes_genus)
rownames(microbes_genus) <- microbes_genus$toAgg
microbes_genus <- microbes_genus[,2:length(microbes_genus)]
microbes_genus <- as.data.frame(t(microbes_genus))
expcn <- colnames(microbes_genus); expcn<-gsub("g__","Genus__",expcn,fixed = T);colnames(microbes_genus)<-expcn

#family
microbes_toagg <- as.data.frame(t(microbes))
stopifnot(all(rownames(microbes_toagg)==rownames(cn)))
microbes_toagg <- cbind(cn$family,microbes_toagg); cnn <- colnames(microbes_toagg); cnn[1] <- "toAgg"; colnames(microbes_toagg) <- cnn
microbes_toagg <- microbes_toagg[microbes_toagg$toAgg!="f__",]

microbes_family <- microbes_toagg %>% 
  group_by(toAgg) %>% 
  summarise(across(everything(), ~sum(.x,na.rm = TRUE)))
microbes_family <- as.data.frame(microbes_family)
rownames(microbes_family) <- microbes_family$toAgg
microbes_family <- microbes_family[,2:length(microbes_family)]
microbes_family <- as.data.frame(t(microbes_family))
expcn <- colnames(microbes_family); expcn<-gsub("f__","Family__",expcn,fixed = T);colnames(microbes_family)<-expcn

##############################################################################################
# convert to relative abundance
microbes_species <- t(microbes)
microbes_species <- scale(microbes_species,center = F,scale = colSums(microbes_species,na.rm = T))
microbes_species <- as.data.frame(t(microbes_species))

microbes_genus <- t(microbes_genus)
microbes_genus <- scale(microbes_genus,center = F,scale = colSums(microbes_genus,na.rm = T))
microbes_genus <- as.data.frame(t(microbes_genus))

microbes_family <- t(microbes_family)
microbes_family <- scale(microbes_family,center = F,scale = colSums(microbes_family))
microbes_family <- as.data.frame(t(microbes_family))

microbes_combined_test <- cbind(microbes_species,microbes_genus,microbes_family)
stopifnot(all(rownames(microbes_combined_test)==rownames(extendedmodel_test)))

##############################################################################################

#load only significant bacteria from association analysis
  
#read results GOLD file
GOLDres <- read.csv(file=paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Bacteria_",TRAIN_SPLIT,"_pval_v4.txt"),sep="\t")
#drop all rnkmrm
GOLDres <- GOLDres[,!grepl("FDR",colnames(GOLDres),fixed = T)]
# c2 <- colnames(GOLDres)[grepl("C2",colnames(GOLDres),fixed = T)]
# c2 <- c2[!grepl("_at30_",c2,fixed = T)]
# GOLDres <- GOLDres[,!colnames(GOLDres) %in% c2]

#filter GOLD by bacteria with at least any of the OSA measurements significant (0.05) which are also in RTOSA
GOLDres_rn <- c()
GOLDres_rnList <- list()

c0 <- colnames(GOLDres)[grepl("AHI_permPval",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[1]] <- unique(c(rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, unique(c(rownames(a)))))

c0 <- colnames(GOLDres)[grepl("MinSpO2_permPval",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[2]] <- unique(c(rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))

c0 <- colnames(GOLDres)[grepl("AvgSpO2_permPval",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[3]] <- unique(c(rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))

c0 <- colnames(GOLDres)[grepl("SpO290_permPval",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[4]] <- unique(c(rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))


c0 <- colnames(GOLDres)[grepl("AHI_C2_at5_perm",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# GOLDres_rnList[[5]] <- unique(c(rownames(a)))
GOLDres_rnList[[5]] <- unique(c(rownames(a),GOLDres_rnList[[1]]))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))

c0 <- colnames(GOLDres)[grepl("AHI_C2_at15_perm",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# GOLDres_rnList[[6]] <- unique(c(rownames(a)))
GOLDres_rnList[[6]] <- unique(c(rownames(a),GOLDres_rnList[[1]]))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))

c0 <- colnames(GOLDres)[grepl("AHI_C2_at30_perm",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# GOLDres_rnList[[7]] <- unique(c(rownames(a)))
GOLDres_rnList[[7]] <- unique(c(rownames(a),GOLDres_rnList[[1]]))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))

c0 <- colnames(GOLDres)[grepl("AHI_C2_at5_15_perm",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# GOLDres_rnList[[8]] <- unique(c(rownames(a)))
GOLDres_rnList[[8]] <- unique(c(rownames(a),GOLDres_rnList[[1]]))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))

c0 <- colnames(GOLDres)[grepl("AHI_C2_at5_30_perm",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# GOLDres_rnList[[9]] <- unique(c(rownames(a)))
GOLDres_rnList[[9]] <- unique(c(rownames(a),GOLDres_rnList[[1]]))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))

microbes_combined <- microbes_combined[,GOLDres_rn]
microbes_combined_test <- microbes_combined_test[,GOLDres_rn]


stopifnot(all(rownames(microbes_combined)==rownames(extendedmodel)))
stopifnot(all(rownames(microbes_combined_test)==rownames(extendedmodel_test)))

###############################################################################################
##############################################################################################

functions_combined <- cbind(pathways,pathwaysKEGGModule,pathwaysKO)
functions_combined_test <- cbind(pathways_test,pathwaysKEGGModule_test,pathwaysKO_test)

#load only significant function from association analysis

#read results GOLD file
FuncGOLDres <- read.csv(file=paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Function_",TRAIN_SPLIT,"_pval_v4.txt"),sep="\t")
#drop all rnkmrm
FuncGOLDres <- FuncGOLDres[,!grepl("FDR",colnames(FuncGOLDres),fixed = T)]
# c2 <- colnames(FuncGOLDres)[grepl("C2",colnames(FuncGOLDres),fixed = T)]
# c2 <- c2[!grepl("_at30_",c2,fixed = T)]
# FuncGOLDres <- FuncGOLDres[,!colnames(FuncGOLDres) %in% c2]

#filter GOLD by bacteria with at least any of the OSA measurements significant (0.05) which are also in RTOSA
FuncGOLDres_rn <- c()
FuncGOLDres_rnList <- list()

c0 <- colnames(FuncGOLDres)[grepl("AHI_permPval",colnames(FuncGOLDres),fixed = T)]
a <- FuncGOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
FuncGOLDres_rnList[[1]] <- unique(c(rownames(a)))
FuncGOLDres_rn <- unique(c(FuncGOLDres_rn, unique(c(rownames(a)))))

c0 <- colnames(FuncGOLDres)[grepl("MinSpO2_permPval",colnames(FuncGOLDres),fixed = T)]
a <- FuncGOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
FuncGOLDres_rnList[[2]] <- unique(c(rownames(a)))
FuncGOLDres_rn <- unique(c(FuncGOLDres_rn, rownames(a)))

c0 <- colnames(FuncGOLDres)[grepl("AvgSpO2_permPval",colnames(FuncGOLDres),fixed = T)]
a <- FuncGOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
FuncGOLDres_rnList[[3]] <- unique(c(rownames(a)))
FuncGOLDres_rn <- unique(c(FuncGOLDres_rn, rownames(a)))

c0 <- colnames(FuncGOLDres)[grepl("SpO290_permPval",colnames(FuncGOLDres),fixed = T)]
a <- FuncGOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
FuncGOLDres_rnList[[4]] <- unique(c(rownames(a)))
FuncGOLDres_rn <- unique(c(FuncGOLDres_rn, rownames(a)))


c0 <- colnames(FuncGOLDres)[grepl("AHI_C2_at5_perm",colnames(FuncGOLDres),fixed = T)]
a <- FuncGOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# FuncGOLDres_rnList[[5]] <- unique(c(rownames(a)))
FuncGOLDres_rnList[[5]] <- unique(c(rownames(a),FuncGOLDres_rnList[[1]]))
FuncGOLDres_rn <- unique(c(FuncGOLDres_rn, rownames(a)))

c0 <- colnames(FuncGOLDres)[grepl("AHI_C2_at15_perm",colnames(FuncGOLDres),fixed = T)]
a <- FuncGOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# FuncGOLDres_rnList[[6]] <- unique(c(rownames(a)))
FuncGOLDres_rnList[[6]] <- unique(c(rownames(a),FuncGOLDres_rnList[[1]]))
FuncGOLDres_rn <- unique(c(FuncGOLDres_rn, rownames(a)))

c0 <- colnames(FuncGOLDres)[grepl("AHI_C2_at30_perm",colnames(FuncGOLDres),fixed = T)]
a <- FuncGOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# FuncGOLDres_rnList[[7]] <- unique(c(rownames(a)))
FuncGOLDres_rnList[[7]] <- unique(c(rownames(a),FuncGOLDres_rnList[[1]]))
FuncGOLDres_rn <- unique(c(FuncGOLDres_rn, rownames(a)))

c0 <- colnames(FuncGOLDres)[grepl("AHI_C2_at5_15_perm",colnames(FuncGOLDres),fixed = T)]
a <- FuncGOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# FuncGOLDres_rnList[[8]] <- unique(c(rownames(a)))
FuncGOLDres_rnList[[8]] <- unique(c(rownames(a),FuncGOLDres_rnList[[1]]))
FuncGOLDres_rn <- unique(c(FuncGOLDres_rn, rownames(a)))

c0 <- colnames(FuncGOLDres)[grepl("AHI_C2_at5_30_perm",colnames(FuncGOLDres),fixed = T)]
a <- FuncGOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
# FuncGOLDres_rnList[[9]] <- unique(c(rownames(a)))
FuncGOLDres_rnList[[9]] <- unique(c(rownames(a),FuncGOLDres_rnList[[1]]))
FuncGOLDres_rn <- unique(c(FuncGOLDres_rn, rownames(a)))

stopifnot(all(rownames(pathways)==rownames(pathwaysKEGGModule)))
stopifnot(all(rownames(pathways)==rownames(pathwaysKO)))

functions_combined <- functions_combined[,FuncGOLDres_rn]
functions_combined_test <- functions_combined_test[,FuncGOLDres_rn]


# real_microbe_names <- colnames(microbes_combined)
# colnames(microbes_combined) <- paste0("bact",1:length(colnames(microbes_combined)))

bidx <- which(rownames(functions_combined)!=rownames(extendedmodel))
stopifnot(length(bidx) < 2)
if(length(bidx)>0){
  functions_combined <- functions_combined[-bidx,]
  extendedmodel <- extendedmodel[-bidx,]
  microbes_combined <- microbes_combined[-bidx,]
}
stopifnot(all(rownames(functions_combined)==rownames(extendedmodel)))


bidx <- which(rownames(functions_combined_test)!=rownames(extendedmodel_test))
stopifnot(length(bidx) < 2)
if(length(bidx)>0){
  functions_combined_test <- functions_combined_test[-bidx,]
  extendedmodel_test <- extendedmodel_test[-bidx,]
}
stopifnot(all(rownames(functions_combined_test)==rownames(extendedmodel_test)))

bidx <- which(!rownames(microbes_combined_test) %in% rownames(extendedmodel_test))
microbes_combined_test <- microbes_combined_test[-bidx,]
stopifnot(all(rownames(microbes_combined_test)==rownames(extendedmodel_test)))

###############################################################################################

outcomes <- c("AHI","MinSpO2","AvgSpO2","SpO290")
outcomesFactors <- c("AHI_C2_at15","AHI_C2_at30","AHI_C2_at5_15","AHI_C2_at5_30")

curmodel <- extendedmodel
curmodel_test <- extendedmodel_test

# curmodel <- mainmodel
#define categories for logistic regression
curmodel$AHI_C2_at5 <- ifelse(curmodel$AHI<5,"no_OSA","OSA")
curmodel$AHI_C2_at15 <- ifelse(curmodel$AHI<15,"no_OSA","OSA")
curmodel$AHI_C2_at30 <- ifelse(curmodel$AHI<30,"no_OSA","OSA")
curmodel$AHI_C2_at5_15 <- ifelse(curmodel$AHI<5,"no_OSA",ifelse(curmodel$AHI>15,"OSA",NA))
curmodel$AHI_C2_at5_30 <- ifelse(curmodel$AHI<5,"no_OSA",ifelse(curmodel$AHI>30,"OSA",NA))

curmodel$AHI_C2_at5 <- ifelse(is.na(curmodel$AHI_C2_at5),"no_OSA",curmodel$AHI_C2_at5)
curmodel$AHI_C2_at15 <- ifelse(is.na(curmodel$AHI_C2_at15),"no_OSA",curmodel$AHI_C2_at15)
curmodel$AHI_C2_at30 <- ifelse(is.na(curmodel$AHI_C2_at30),"no_OSA",curmodel$AHI_C2_at30)
curmodel$AHI_C2_at5_15 <- ifelse(is.na(curmodel$AHI_C2_at5_15),"no_OSA",curmodel$AHI_C2_at5_15)
curmodel$AHI_C2_at5_30 <- ifelse(is.na(curmodel$AHI_C2_at5_30),"no_OSA",curmodel$AHI_C2_at5_30)

curmodel_test$AHI_C2_at5 <- ifelse(curmodel_test$AHI<5,"no_OSA","OSA")
curmodel_test$AHI_C2_at15 <- ifelse(curmodel_test$AHI<15,"no_OSA","OSA")
curmodel_test$AHI_C2_at30 <- ifelse(curmodel_test$AHI<30,"no_OSA","OSA")
curmodel_test$AHI_C2_at5_15 <- ifelse(curmodel_test$AHI<5,"no_OSA",ifelse(curmodel_test$AHI>15,"OSA",NA))
curmodel_test$AHI_C2_at5_30 <- ifelse(curmodel_test$AHI<5,"no_OSA",ifelse(curmodel_test$AHI>30,"OSA",NA))

curmodel_test$AHI_C2_at5 <- ifelse(is.na(curmodel_test$AHI_C2_at5),"no_OSA",curmodel_test$AHI_C2_at5)
curmodel_test$AHI_C2_at15 <- ifelse(is.na(curmodel_test$AHI_C2_at15),"no_OSA",curmodel_test$AHI_C2_at15)
curmodel_test$AHI_C2_at30 <- ifelse(is.na(curmodel_test$AHI_C2_at30),"no_OSA",curmodel_test$AHI_C2_at30)
curmodel_test$AHI_C2_at5_15 <- ifelse(is.na(curmodel_test$AHI_C2_at5_15),"no_OSA",curmodel_test$AHI_C2_at5_15)
curmodel_test$AHI_C2_at5_30 <- ifelse(is.na(curmodel_test$AHI_C2_at5_30),"no_OSA",curmodel_test$AHI_C2_at5_30)

#define models, covariates
# covariates <- c("AGE","Sex","CENTER","WEIGHT_NORM_OVERALL_V2")
covariates <- c("AGE","Sex","CENTER","WEIGHT_NORM_OVERALL_V2",
                "EDUCATION_C3_V2",
                "INCOME_C5_V2",                               
                "CURRENT_SMOKER_V2","ALCOHOL_INTAKE_V2",
                "HYPERTENSION_C4_V2",
                "DIABETES4_V2",
                "BMI_V2","WAIST_HIP_V2","mediter_diet_score")

####################################### bacterial ######################################################
names(GOLDres_rnList) <- c(outcomes, outcomesFactors)

bactres <- matrix(nrow = length(outcomesFactors),ncol = 7)
rownames(bactres) <- outcomesFactors; colnames(bactres) <- c("Covariates LM","Covariates LASSO",
                                                             "Microbes LM","Microbes LASSO",
                                                             "Covariates & Microbes LM", "Covariates & Microbes LASSO",
                                                             "All Covariates and LASSO Microbes (LM)")
bactSelected <- bactres
bactModels <- bactres
bactModels <- as.data.frame(bactModels)

data_ctrl <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary,
                          classProbs = TRUE)
set.seed(42)
for(outcome in outcomesFactors){ #outcome <- outcomesFactors[2]
  print(paste("Weighted classification : ",outcome))

  ##################### Covariates only - LASSO ###################################
  tdata <- cbind(curmodel[,c(outcomesFactors,covariates)])
  tdata[[outcome]] <- as.factor(tdata[[outcome]])
  model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(covariates),collapse="+")))
  real_microbe_names <- colnames(tdata)
  real_microbe_names <- as.data.frame(real_microbe_names)
  real_microbe_names$Encoding <- colnames(tdata)
  rownames(real_microbe_names) <- real_microbe_names$Encoding
  
  
  model_weights <- ifelse(tdata[[outcome]] == "no_OSA",
                          (1/table(tdata[[outcome]])[1]) * 0.5,
                          (1/table(tdata[[outcome]])[2]) * 0.5)
  
  model_caretLASSO <- train(model.formula,   # model to fit
                       data = tdata,                        
                       trControl = data_ctrl,              # folds
                       method = "glmnet",                      # specifying regression model
                       tuneGrid=expand.grid(
                         .alpha=1,
                         .lambda=seq(0, 10, by = 0.1)),
                       family = "binomial",
                       weights = model_weights,
                       metric = "ROC",
                       na.action = na.omit)                # pass missing data to model - some models will handle this
  bestcoeff <- data.frame(as.matrix(coef(model_caretLASSO$finalModel, model_caretLASSO$bestTune$lambda)))
  nonzerofunc <- rownames(bestcoeff)[which(bestcoeff!=0)]; nonzerofunc <- nonzerofunc[!grepl("(Intercept)",nonzerofunc,fixed = T)]

  #test on the validation set  
  ddata <- cbind(curmodel_test[,c(outcomesFactors,covariates)])
  ddata[[outcome]] <- as.factor(ddata[[outcome]])
  ddata <- na.omit(ddata)
  model_pred <- predict(model_caretLASSO, ddata,type = "prob", verbose=F)
  result.roc <- roc(ddata[[outcome]], model_pred$OSA,
                    smoothed = T,auc = T,
                    ci=TRUE, ci.alpha=0.95, quiet = T)

  nz <- c()
  for(n in nonzerofunc){ # n = nonzerofunc[2]
    for(m in real_microbe_names$Encoding){
      if(grepl(m,n,fixed = T)){nz <- c(nz,m)}
    }
  }
  nz <- unique(nz)
  bactSelected[outcome,"Covariates LASSO"] <- paste0(real_microbe_names[nz,"real_microbe_names"],collapse = ",")
  bactres[outcome,"Covariates LASSO"] <- result.roc$auc
  bactModels[outcome,"Covariates LASSO"] <- list(list(result.roc))

  print(paste0("Covariates LASSO AUC: ",round(result.roc$auc,2)," Selected ",paste0(real_microbe_names[nz,"real_microbe_names"],collapse = ",")))

  ##################### Covariates only - LM ###################################
  model_caretLM <- train(model.formula,   # model to fit
                            data = tdata,                        
                            trControl = data_ctrl,              # folds
                            method = "glm",                      # specifying regression model
                            family = "binomial",
                            weights = model_weights,
                            metric = "ROC",
                            na.action = na.omit)                # pass missing data to model - some models will handle this

  #test on the validation set  
  model_predLM <- predict(model_caretLM, ddata,type = "prob", verbose=F)
  result.rocLM <- roc(ddata[[outcome]], model_predLM$OSA, quiet = T) # Draw ROC curve.
  print(paste0("Covariates LM AUC: ",round(result.rocLM$auc,2)))
  bactres[outcome,"Covariates LM"] <- result.rocLM$auc
  bactSelected[outcome,"Covariates LM"] <- paste0(covariates,collapse = ",")
  
  
  
  ################################### Microbes only LASSO #####################################3
  tdata <- microbes_combined[,GOLDres_rnList[[outcome]]]
  tdata <- tdata[, colSums(tdata != 0) > 0]
  real_microbe_names <- colnames(tdata)
  colnames(tdata) <- paste0("bact",1:length(real_microbe_names))
  real_microbe_names <- as.data.frame(real_microbe_names)
  real_microbe_names$Encoding <- colnames(tdata)
  rownames(real_microbe_names) <- real_microbe_names$Encoding
  
  tdata <- cbind(curmodel[,c(outcome)],tdata)
  colnames(tdata) <- c(outcome, colnames(tdata)[2:length(tdata)])
  tdata[[outcome]] <- as.factor(tdata[[outcome]])
  model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(paste0("bact",1:length(real_microbe_names))),collapse="+")))
  
  model_caretLASSO <- train(model.formula,   # model to fit
                       data = tdata,                        
                       trControl = data_ctrl,              # folds
                       method = "glmnet",                      # specifying regression model
                       tuneGrid=expand.grid(
                         .alpha=1,
                         .lambda=seq(0, 10, by = 0.1)),
                       family = "binomial",
                       weights = model_weights,
                       metric = "ROC",
                       na.action = na.omit,
                       verbose=F)                # pass missing data to model - some models will handle this
  bestcoeff <- data.frame(as.matrix(coef(model_caretLASSO$finalModel, model_caretLASSO$bestTune$lambda)))
  nonzerofunc <- rownames(bestcoeff)[which(bestcoeff!=0)]; nonzerofunc <- nonzerofunc[!grepl("(Intercept)",nonzerofunc,fixed = T)]

  #test on the validation set  
  ddata <- microbes_combined_test[,GOLDres_rnList[[outcome]]]
  real_microbe_names_test <- colnames(ddata)
  colnames(ddata) <- paste0("bact",1:length(real_microbe_names_test))
  
  ddata <- cbind(curmodel_test[,c(outcome)],ddata)
  colnames(ddata) <- c(outcome, colnames(ddata)[2:length(ddata)])
  ddata[[outcome]] <- as.factor(ddata[[outcome]])
  
  ddata <- na.omit(ddata)
  model_pred <- predict(model_caretLASSO, ddata,type = "prob", verbose=F)
  result.roc <- roc(ddata[[outcome]], model_pred$OSA,
                    smoothed = T,auc = T,
                    ci=TRUE, ci.alpha=0.95, quiet = T)
  print(paste0("Microbes LASSO AUC: ",round(result.roc$auc,2)," Selected ",paste0(real_microbe_names[nonzerofunc,"real_microbe_names"],collapse = ",")))
  bactres[outcome,"Microbes LASSO"] <- result.roc$auc
  bactSelected[outcome,"Microbes LASSO"] <- paste0(real_microbe_names[nonzerofunc,"real_microbe_names"],collapse = ",")
  bactModels[outcome,"Microbes LASSO"] <- list(list(result.roc))
  

  ##################### Microbes only - LM ###################################
  model_caretLM <- train(model.formula,   # model to fit
                         data = tdata,                        
                         trControl = data_ctrl,              # folds
                         method = "glm",                      # specifying regression model
                         family = "binomial",
                         weights = model_weights,
                         metric = "ROC",
                         na.action = na.omit)                # pass missing data to model - some models will handle this
  
  #test on the validation set  
  model_predLM <- predict(model_caretLM, ddata,type = "prob", verbose=F)
  result.rocLM <- roc(ddata[[outcome]], model_predLM$OSA, quiet = T) # Draw ROC curve.
  print(paste0("Microbes LM AUC: ",round(result.rocLM$auc,2)))
  bactres[outcome,"Microbes LM"] <- result.rocLM$auc
  bactSelected[outcome,"Microbes LM"] <- paste0(real_microbe_names[,"real_microbe_names"],collapse = ",")
  
  
  ###########################Covariates and Microbes LASSO#########################################################
  tdata <- microbes_combined[,GOLDres_rnList[[outcome]]]
  tdata <- tdata[, colSums(tdata != 0) > 0]
  real_microbe_names <- colnames(tdata)
  colnames(tdata) <- paste0("bact",1:length(real_microbe_names))

  tdata <- cbind(curmodel[,c(outcomesFactors,covariates)],tdata)
  tdata[[outcome]] <- as.factor(tdata[[outcome]])
  model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(covariates,paste0("bact",1:length(real_microbe_names))),collapse="+")))

  real_microbe_names <- c(outcomesFactors,covariates,real_microbe_names)
  real_microbe_names <- as.data.frame(real_microbe_names)
  real_microbe_names$Encoding <- colnames(tdata)
  rownames(real_microbe_names) <- real_microbe_names$Encoding
  
  model_caretLASSO <- train(model.formula,   # model to fit
                       data = tdata,                        
                       trControl = data_ctrl,              # folds
                       method = "glmnet",                      # specifying regression model
                       tuneGrid=expand.grid(
                         .alpha=1,
                         .lambda=seq(0, 10, by = 0.1)),
                       family = "binomial",
                       weights = model_weights,
                       metric = "ROC",
                       na.action = na.omit)                # pass missing data to model - some models will handle this
  bestcoeff <- data.frame(as.matrix(coef(model_caretLASSO$finalModel, model_caretLASSO$bestTune$lambda)))
  nonzerofunc <- rownames(bestcoeff)[which(bestcoeff!=0)]; nonzerofunc <- nonzerofunc[!grepl("(Intercept)",nonzerofunc,fixed = T)]
  
  #test on the validation set  
  ddata <- microbes_combined_test[,GOLDres_rnList[[outcome]]]
  real_microbe_names_test <- colnames(ddata)
  colnames(ddata) <- paste0("bact",1:length(real_microbe_names_test))

  ddata <- cbind(curmodel_test[,c(outcome,covariates)],ddata)
  colnames(ddata) <- c(outcome, colnames(ddata)[2:length(ddata)])
  ddata[[outcome]] <- as.factor(ddata[[outcome]])
  
  ddata <- na.omit(ddata)
  model_pred <- predict(model_caretLASSO, ddata,type = "prob")
  result.roc <- roc(ddata[[outcome]], model_pred$OSA,
                    smoothed = T,auc = T,
                    ci=TRUE, ci.alpha=0.95, quiet = T) # Draw ROC curve.
  
  nz <- c()
  for(n in nonzerofunc){ # n = nonzerofunc[2]
    for(m in real_microbe_names$Encoding){
      if(grepl(m,n,fixed = T)){nz <- c(nz,m)}
    }
  }
  nz <- unique(nz)

  print(paste0("Covariates & Microbes LASSO AUC: ",round(result.roc$auc,2)," Selected ",paste0(real_microbe_names[nz,"real_microbe_names"],collapse = ",")))
  bactres[outcome,"Covariates & Microbes LASSO"] <- result.roc$auc
  bactSelected[outcome,"Covariates & Microbes LASSO"] <- paste0(real_microbe_names[nz,"real_microbe_names"],collapse = ",")
  bactModels[outcome,"Covariates & Microbes LASSO"] <- list(list(result.roc))
  
  ##################### Covariates & Microbes - LM ###################################
  model_caretLM <- train(model.formula,   # model to fit
                         data = tdata,                        
                         trControl = data_ctrl,              # folds
                         method = "glm",                      # specifying regression model
                         family = "binomial",
                         weights = model_weights,
                         metric = "ROC",
                         na.action = na.omit)                # pass missing data to model - some models will handle this
  
  #test on the validation set  
  model_predLM <- predict(model_caretLM, ddata,type = "prob", verbose=F)
  result.rocLM <- roc(ddata[[outcome]], model_predLM$OSA, quiet = T) # Draw ROC curve.
  print(paste0("Covariates & Microbes LM AUC: ",round(result.rocLM$auc,2)))
  bactres[outcome,"Covariates & Microbes LM"] <- result.rocLM$auc
  bactSelected[outcome,"Covariates & Microbes LM"] <- paste0(real_microbe_names[,"real_microbe_names"],collapse = ",")
  
  ###########################All Covariates and LASSO Microbes (LM) #########################################################
  if(nchar(bactSelected[outcome,"Microbes LASSO"])>0){
    tdata <- as.data.frame(microbes_combined[,strsplit(bactSelected[outcome,"Microbes LASSO"],",")[[1]]])
    real_microbe_names <- colnames(tdata)
    colnames(tdata) <- paste0("bact",1:length(real_microbe_names))
  
    tdata <- cbind(curmodel[,c(outcomesFactors,covariates)],tdata)
    tdata[[outcome]] <- as.factor(tdata[[outcome]])
    model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(covariates,paste0("bact",1:length(real_microbe_names))),collapse="+")))
    
    real_microbe_names <- c(outcomesFactors,covariates,real_microbe_names)
    real_microbe_names <- as.data.frame(real_microbe_names)
    real_microbe_names$Encoding <- colnames(tdata)
    rownames(real_microbe_names) <- real_microbe_names$Encoding
  }else{
    print(paste0("All Covariates and LASSO Microbes (LM): NA"))
    bactres[outcome,"All Covariates and LASSO Microbes (LM)"] <- NA
    bactSelected[outcome,"All Covariates and LASSO Microbes (LM)"] <- ""
    next
  }
  
  
  model_caretLM <- train(model.formula,   # model to fit
                            data = tdata,                        
                            trControl = data_ctrl,              # folds
                            method = "glm",                      # specifying regression model
                            family = "binomial",
                            weights = model_weights,
                            metric = "ROC",
                            na.action = na.omit)                # pass missing data to model - some models will handle this

  #test on the validation set  
  ddata <- as.data.frame(microbes_combined_test[,strsplit(bactSelected[outcome,"Microbes LASSO"],",")[[1]]])
  real_microbe_names_test <- colnames(ddata)
  colnames(ddata) <- paste0("bact",1:length(real_microbe_names_test))
  
  ddata <- cbind(curmodel_test[,c(outcome,covariates)],ddata)
  colnames(ddata) <- c(outcome, colnames(ddata)[2:length(ddata)])
  ddata[[outcome]] <- as.factor(ddata[[outcome]])
  
  ddata <- na.omit(ddata)
  model_pred <- predict(model_caretLM, ddata,type = "prob")
  result.roc <- roc(ddata[[outcome]], model_pred$OSA, quiet = T) # Draw ROC curve.
  print(paste0("All Covariates and LASSO Microbes (LM): ",round(result.roc$auc,2)))
  bactres[outcome,"All Covariates and LASSO Microbes (LM)"] <- result.roc$auc
  bactSelected[outcome,"All Covariates and LASSO Microbes (LM)"] <- paste0(real_microbe_names[,"real_microbe_names"],collapse = ",")
  
}
print(bactres)

if(ressave){
  write.table(bactres,
              file=paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Bacteria_GOLD_BRS_AUC_",TRAIN_SPLIT,ressuf,".txt"),
              sep="\t")
  write.table(bactSelected,
              file=paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Bacteria_GOLD_BRS_LASSOSelected_",TRAIN_SPLIT,ressuf,".txt"),
              sep="\t")
  write.table(apply(bactSelected,c(1,2),function(x){length(strsplit(x,",",fixed = T)[[1]])}),
              file=paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Bacteria_GOLD_BRS_LASSOSelectedNum_",TRAIN_SPLIT,ressuf,".txt"),
              sep="\t")
}

########################################## functional ###################################################
names(FuncGOLDres_rnList) <- c(outcomes, outcomesFactors)


funcres <- matrix(nrow = length(outcomesFactors),ncol = 7)
rownames(funcres) <- outcomesFactors; colnames(funcres) <- c("Covariates LM","Covariates LASSO",
                                                             "Functions LM","Functions LASSO",
                                                             "Covariates & Functions LM", "Covariates & Functions LASSO",
                                                             "All Covariates and LASSO Functions (LM)")

funcSelected <- funcres
funcModels <- funcres; funcModels <- as.data.frame(funcModels)
data_ctrl <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary,
                          classProbs = TRUE)
set.seed(42)
for(outcome in outcomesFactors){ #outcome <- outcomesFactors[2]
  print(paste("Weighted classification : ",outcome))
  
  ##################### Covariates only - LASSO ###################################
  tdata <- cbind(curmodel[,c(outcomesFactors,covariates)])
  tdata[[outcome]] <- as.factor(tdata[[outcome]])
  model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(covariates),collapse="+")))
  real_microbe_names <- colnames(tdata)
  real_microbe_names <- as.data.frame(real_microbe_names)
  real_microbe_names$Encoding <- colnames(tdata)
  rownames(real_microbe_names) <- real_microbe_names$Encoding
  
  
  model_weights <- ifelse(tdata[[outcome]] == "no_OSA",
                          (1/table(tdata[[outcome]])[1]) * 0.5,
                          (1/table(tdata[[outcome]])[2]) * 0.5)
  
  model_caretLASSO <- train(model.formula,   # model to fit
                            data = tdata,                        
                            trControl = data_ctrl,              # folds
                            method = "glmnet",                      # specifying regression model
                            tuneGrid=expand.grid(
                              .alpha=1,
                              .lambda=seq(0, 10, by = 0.1)),
                            family = "binomial",
                            weights = model_weights,
                            metric = "ROC",
                            na.action = na.omit)                # pass missing data to model - some models will handle this
  bestcoeff <- data.frame(as.matrix(coef(model_caretLASSO$finalModel, model_caretLASSO$bestTune$lambda)))
  nonzerofunc <- rownames(bestcoeff)[which(bestcoeff!=0)]; nonzerofunc <- nonzerofunc[!grepl("(Intercept)",nonzerofunc,fixed = T)]
  
  #test on the validation set  
  ddata <- cbind(curmodel_test[,c(outcomesFactors,covariates)])
  ddata[[outcome]] <- as.factor(ddata[[outcome]])
  ddata <- na.omit(ddata)
  model_pred <- predict(model_caretLASSO, ddata,type = "prob", verbose=F)
  result.roc <- roc(ddata[[outcome]], model_pred$OSA) # Draw ROC curve.
  
  nz <- c()
  for(n in nonzerofunc){ # n = nonzerofunc[2]
    for(m in real_microbe_names$Encoding){
      if(grepl(m,n,fixed = T)){nz <- c(nz,m)}
    }
  }
  nz <- unique(nz)
  funcSelected[outcome,"Covariates LASSO"] <- paste0(real_microbe_names[nz,"real_microbe_names"],collapse = "//")
  funcres[outcome,"Covariates LASSO"] <- result.roc$auc
  print(paste0("Covariates LASSO AUC: ",round(result.roc$auc,2)," Selected ",paste0(real_microbe_names[nz,"real_microbe_names"],collapse = "//")))
  funcModels[outcome,"Covariates LASSO"] <- list(list(result.roc))
  
  ##################### Covariates only - LM ###################################
  model_caretLM <- train(model.formula,   # model to fit
                         data = tdata,                        
                         trControl = data_ctrl,              # folds
                         method = "glm",                      # specifying regression model
                         family = "binomial",
                         weights = model_weights,
                         metric = "ROC",
                         na.action = na.omit)                # pass missing data to model - some models will handle this
  
  #test on the validation set  
  model_predLM <- predict(model_caretLM, ddata,type = "prob", verbose=F)
  result.rocLM <- roc(ddata[[outcome]], model_predLM$OSA) # Draw ROC curve.
  print(paste0("Covariates LM AUC: ",round(result.rocLM$auc,2)))
  funcres[outcome,"Covariates LM"] <- result.rocLM$auc
  funcSelected[outcome,"Covariates LM"] <- paste0(covariates,collapse = "//")
  
  
  
  ################################### Functions only LASSO #####################################3
  tdata <- functions_combined[,FuncGOLDres_rnList[[outcome]]]
  tdata <- tdata[, colSums(tdata != 0) > 0]
  real_microbe_names <- colnames(tdata)
  colnames(tdata) <- paste0("func",1:length(real_microbe_names))
  real_microbe_names <- as.data.frame(real_microbe_names)
  real_microbe_names$Encoding <- colnames(tdata)
  rownames(real_microbe_names) <- real_microbe_names$Encoding
  
  tdata <- cbind(curmodel[,c(outcome)],tdata)
  colnames(tdata) <- c(outcome, colnames(tdata)[2:length(tdata)])
  tdata[[outcome]] <- as.factor(tdata[[outcome]])
  model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(paste0("func",1:length(real_microbe_names))),collapse="+")))
  
  model_caretLASSO <- train(model.formula,   # model to fit
                            data = tdata,                        
                            trControl = data_ctrl,              # folds
                            method = "glmnet",                      # specifying regression model
                            tuneGrid=expand.grid(
                              .alpha=1,
                              .lambda=seq(0, 10, by = 0.1)),
                            family = "binomial",
                            weights = model_weights,
                            metric = "ROC",
                            na.action = na.omit,
                            verbose=F)                # pass missing data to model - some models will handle this
  bestcoeff <- data.frame(as.matrix(coef(model_caretLASSO$finalModel, model_caretLASSO$bestTune$lambda)))
  nonzerofunc <- rownames(bestcoeff)[which(bestcoeff!=0)]; nonzerofunc <- nonzerofunc[!grepl("(Intercept)",nonzerofunc,fixed = T)]
  
  #test on the validation set  
  ddata <- functions_combined_test[,FuncGOLDres_rnList[[outcome]]]
  real_microbe_names_test <- colnames(ddata)
  colnames(ddata) <- paste0("func",1:length(real_microbe_names_test))
  
  ddata <- cbind(curmodel_test[,c(outcome)],ddata)
  colnames(ddata) <- c(outcome, colnames(ddata)[2:length(ddata)])
  ddata[[outcome]] <- as.factor(ddata[[outcome]])
  
  ddata <- na.omit(ddata)
  model_pred <- predict(model_caretLASSO, ddata,type = "prob", verbose=F)
  result.roc <- roc(ddata[[outcome]], model_pred$OSA) # Draw ROC curve.
  print(paste0("Functions LASSO AUC: ",round(result.roc$auc,2)," Selected ",paste0(real_microbe_names[nonzerofunc,"real_microbe_names"],collapse = "//")))
  funcres[outcome,"Functions LASSO"] <- result.roc$auc
  funcSelected[outcome,"Functions LASSO"] <- paste0(real_microbe_names[nonzerofunc,"real_microbe_names"],collapse = "//")
  funcModels[outcome,"Functions LASSO"] <- list(list(result.roc))
  
  ##################### Functions only - LM ###################################
  model_caretLM <- train(model.formula,   # model to fit
                         data = tdata,                        
                         trControl = data_ctrl,              # folds
                         method = "glm",                      # specifying regression model
                         family = "binomial",
                         weights = model_weights,
                         metric = "ROC",
                         na.action = na.omit)                # pass missing data to model - some models will handle this
  
  #test on the validation set  
  model_predLM <- predict(model_caretLM, ddata,type = "prob", verbose=F)
  result.rocLM <- roc(ddata[[outcome]], model_predLM$OSA) # Draw ROC curve.
  print(paste0("Functions LM AUC: ",round(result.rocLM$auc,2)))
  funcres[outcome,"Functions LM"] <- result.rocLM$auc
  funcSelected[outcome,"Functions LM"] <- paste0(real_microbe_names[,"real_microbe_names"],collapse = "//")
  
  
  ###########################Covariates and Functions LASSO#########################################################
  tdata <- functions_combined[,FuncGOLDres_rnList[[outcome]]]
  tdata <- tdata[, colSums(tdata != 0) > 0]
  real_microbe_names <- colnames(tdata)
  colnames(tdata) <- paste0("func",1:length(real_microbe_names))
  
  tdata <- cbind(curmodel[,c(outcomesFactors,covariates)],tdata)
  tdata[[outcome]] <- as.factor(tdata[[outcome]])
  model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(covariates,paste0("func",1:length(real_microbe_names))),collapse="+")))
  
  real_microbe_names <- c(outcomesFactors,covariates,real_microbe_names)
  real_microbe_names <- as.data.frame(real_microbe_names)
  real_microbe_names$Encoding <- colnames(tdata)
  rownames(real_microbe_names) <- real_microbe_names$Encoding
  real_microbe_names <- real_microbe_names[2:dim(real_microbe_names)[1],]
  
  model_caretLASSO <- train(model.formula,   # model to fit
                            data = tdata,                        
                            trControl = data_ctrl,              # folds
                            method = "glmnet",                      # specifying regression model
                            tuneGrid=expand.grid(
                              .alpha=1,
                              .lambda=seq(0, 10, by = 0.1)),
                            family = "binomial",
                            weights = model_weights,
                            metric = "ROC",
                            na.action = na.omit)                # pass missing data to model - some models will handle this
  bestcoeff <- data.frame(as.matrix(coef(model_caretLASSO$finalModel, model_caretLASSO$bestTune$lambda)))
  nonzerofunc <- rownames(bestcoeff)[which(bestcoeff!=0)]; nonzerofunc <- nonzerofunc[!grepl("(Intercept)",nonzerofunc,fixed = T)]
  
  #test on the validation set  
  ddata <- functions_combined_test[,FuncGOLDres_rnList[[outcome]]]
  real_microbe_names_test <- colnames(ddata)
  colnames(ddata) <- paste0("func",1:length(real_microbe_names_test))
  
  ddata <- cbind(curmodel_test[,c(outcome,covariates)],ddata)
  colnames(ddata) <- c(outcome, colnames(ddata)[2:length(ddata)])
  ddata[[outcome]] <- as.factor(ddata[[outcome]])
  
  ddata <- na.omit(ddata)
  model_pred <- predict(model_caretLASSO, ddata,type = "prob")
  result.roc <- roc(ddata[[outcome]], model_pred$OSA) # Draw ROC curve.
  
  nz <- c()
  for(n in nonzerofunc){ # n = nonzerofunc[2]
    for(m in real_microbe_names$Encoding){
      if(grepl(m,n,fixed = T)){nz <- c(nz,m)}
    }
  }
  nz <- unique(nz)
  
  print(paste0("Covariates & Functions LASSO AUC: ",round(result.roc$auc,2)," Selected ",paste0(real_microbe_names[nz,"real_microbe_names"],collapse = "//")))
  funcres[outcome,"Covariates & Functions LASSO"] <- result.roc$auc
  funcSelected[outcome,"Covariates & Functions LASSO"] <- paste0(real_microbe_names[nz,"real_microbe_names"],collapse = "//")
  funcModels[outcome,"Covariates & Functions LASSO"] <- list(list(result.roc))
  
  ##################### Covariates & Functions - LM ###################################
  model_caretLM <- train(model.formula,   # model to fit
                         data = tdata,                        
                         trControl = data_ctrl,              # folds
                         method = "glm",                      # specifying regression model
                         family = "binomial",
                         weights = model_weights,
                         metric = "ROC",
                         na.action = na.omit)                # pass missing data to model - some models will handle this
  
  #test on the validation set  
  model_predLM <- predict(model_caretLM, ddata,type = "prob", verbose=F)
  result.rocLM <- roc(ddata[[outcome]], model_predLM$OSA) # Draw ROC curve.
  print(paste0("Covariates & Functions LM AUC: ",round(result.rocLM$auc,2)))
  funcres[outcome,"Covariates & Functions LM"] <- result.rocLM$auc
  funcSelected[outcome,"Covariates & Functions LM"] <- paste0(real_microbe_names[,"real_microbe_names"],collapse = "//")
  
  ###########################All Covariates and LASSO Functions (LM) #########################################################
  if(nchar(funcSelected[outcome,"Functions LASSO"])>0){
    tdata <- as.data.frame(functions_combined[,strsplit(funcSelected[outcome,"Functions LASSO"],"//")[[1]]])
    real_microbe_names <- colnames(tdata)
    colnames(tdata) <- paste0("func",1:length(real_microbe_names))
    
    tdata <- cbind(curmodel[,c(outcomesFactors,covariates)],tdata)
    tdata[[outcome]] <- as.factor(tdata[[outcome]])
    model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(covariates,paste0("func",1:length(real_microbe_names))),collapse="+")))
    
    real_microbe_names <- c(outcomesFactors,covariates,real_microbe_names)
    real_microbe_names <- as.data.frame(real_microbe_names)
    real_microbe_names$Encoding <- colnames(tdata)
    rownames(real_microbe_names) <- real_microbe_names$Encoding
    real_microbe_names <- real_microbe_names[2:dim(real_microbe_names)[1],]
  }else{
    print(paste0("All Covariates and LASSO Functions (LM): NA"))
    funcres[outcome,"All Covariates and LASSO Functions (LM)"] <- NA
    funcSelected[outcome,"All Covariates and LASSO Functions (LM)"] <- ""
    next
  }
  
  model_caretLM <- train(model.formula,   # model to fit
                         data = tdata,                        
                         trControl = data_ctrl,              # folds
                         method = "glm",                      # specifying regression model
                         family = "binomial",
                         weights = model_weights,
                         metric = "ROC",
                         na.action = na.pass)                # pass missing data to model - some models will handle this
  
  #test on the validation set  
  ddata <- as.data.frame(functions_combined_test[,strsplit(funcSelected[outcome,"Functions LASSO"],"//")[[1]]])
  real_microbe_names_test <- colnames(ddata)
  colnames(ddata) <- paste0("func",1:length(real_microbe_names_test))
  
  ddata <- cbind(curmodel_test[,c(outcome,covariates)],ddata)
  colnames(ddata) <- c(outcome, colnames(ddata)[2:length(ddata)])
  ddata[[outcome]] <- as.factor(ddata[[outcome]])
  
  ddata <- na.omit(ddata)
  model_pred <- predict(model_caretLM, ddata,type = "prob")
  result.roc <- roc(ddata[[outcome]], model_pred$OSA) # Draw ROC curve.
  print(paste0("All Covariates and LASSO Functions (LM): ",round(result.roc$auc,2)))
  funcres[outcome,"All Covariates and LASSO Functions (LM)"] <- result.roc$auc
  funcSelected[outcome,"All Covariates and LASSO Functions (LM)"] <- paste0(real_microbe_names[,"real_microbe_names"],collapse = "//")
  
}
print(funcres)

if(ressave){
  write.table(funcres,
            file=paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Function_GOLD_BRS_AUC_",TRAIN_SPLIT,ressuf,".txt"),
            sep="\t")
  write.table(funcSelected,
            file=paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Function_GOLD_BRS_LASSOSelected_",TRAIN_SPLIT,ressuf,".txt"),
            sep="\t")
  write.table(apply(funcSelected,c(1,2),function(x){length(strsplit(x,"//",fixed = T)[[1]])}),
              file=paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Function_GOLD_BRS_LASSOSelectedNum_",TRAIN_SPLIT,ressuf,".txt"),
              sep="\t")
}
# print(funcSelected)

#############################################################################################################
#plot
outcome <- outcomesFactors[2] # AHI > 30  # outcome2 <- outcomesFactors[1]

roc_GOLD_Covariates_LASSOs <- smooth(bactModels[outcome,"Covariates LASSO"][[1]], method="binormal")
s1 <- ci.sp(roc_GOLD_Covariates_LASSOs,conf.level=0.99,boot.n=10000)
roc_GOLD_Microbes_LASSOs <- smooth(bactModels[outcome,"Microbes LASSO"][[1]], method="binormal",plot=T)
s2 <- ci.sp(roc_GOLD_Microbes_LASSOs,conf.level=0.99)
roc_GOLD_Covariates_Microbes_LASSOs <- smooth(bactModels[outcome,"Covariates & Microbes LASSO"][[1]], method="binormal",plot=T)
s3 <- ci.sp(roc_GOLD_Covariates_Microbes_LASSOs,conf.level=0.99)

rocf_GOLD_Covariates_LASSOs <- smooth(funcModels[outcome,"Covariates LASSO"][[1]], method="binormal")
s7 <- ci.sp(rocf_GOLD_Covariates_LASSOs,conf.level=0.95)
rocf_GOLD_Functions_LASSOs <- smooth(funcModels[outcome,"Functions LASSO"][[1]], method="binormal")
s8 <- ci.sp(rocf_GOLD_Functions_LASSOs,conf.level=0.95)
rocf_GOLD_Covariates_Functions_LASSOs <- smooth(funcModels[outcome,"Covariates & Functions LASSO"][[1]], method="binormal")
s9 <- ci.sp(rocf_GOLD_Covariates_Functions_LASSOs,conf.level=0.95)

if(ressave){
  pdf(paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/BRS_GOLD_RTOSA",ressuf,".pdf"), 
    height=20, width=20)
}

m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.4,0.4,0.2))

par(mar = c(2,2,1,1))
plot.roc(roc_GOLD_Covariates_LASSOs, lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE,col="darkblue", main="GOLD", asp = NA)
plot(s1, type="shape", col="lightblue")
plot.roc(roc_GOLD_Microbes_LASSOs,add=T,col="darkgreen", lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE, print.auc.y=0.4, asp = NA)
plot(s2, type="shape", col="lightgreen")
plot.roc(roc_GOLD_Covariates_Microbes_LASSOs,add=T,col="darkred", lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE, print.auc.y=0.3, asp = NA)

par(mar = c(2,2,1,1))
plot.roc(roc_RTOSA_Covariates_LASSOs, lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE,col="darkblue", main="RTOSA", asp = NA)
plot(s4, type="shape", col="lightblue")
plot.roc(roc_RTOSA_Microbes_LASSO,add=T,col="darkgreen", lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE, print.auc.y=0.4, asp = NA)
plot(s5, type="shape", col="lightgreen")
plot.roc(roc_RTOSA_Covariates_Microbes_LASSOs,add=T,col="darkred", lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE, print.auc.y=0.3, asp = NA)
plot(s6, type="shape", col="pink")

par(mar = c(2,2,1,1))
plot.roc(rocf_GOLD_Covariates_LASSOs, lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE,col="darkblue", main="GOLD", asp = NA)
plot.roc(rocf_GOLD_Functions_LASSOs,add=T,col="darkgreen", lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE, print.auc.y=0.4, asp = NA)
plot.roc(rocf_GOLD_Covariates_Functions_LASSOs,add=T,col="darkred", lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE, print.auc.y=0.3, asp = NA)

par(mar = c(2,2,1,1))
plot.roc(rocf_RTOSA_Covariates_LASSOs, lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE,col="darkblue", main="RTOSA", asp = NA)
plot.roc(rocf_RTOSA_Functions_LASSOs,add=T,col="darkgreen", lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE, print.auc.y=0.4, asp = NA)
plot.roc(rocf_RTOSA_Covariates_Functions_LASSO,add=T,col="darkred", lwd=2, legacy.axes=FALSE, percent=TRUE, print.auc=TRUE, print.auc.y=0.3, asp = NA)

plot(s11, type="shape", col="lightgreen")
plot(s12, type="shape", col="lightblue")

par(mar=c(0,0,1,0))
plot(0.1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("darkblue","darkgreen", "darkred")
legend(x = "top",inset = 0,
       legend=c("Covariates", "Functions","Covariates & Functions"),
       col=plot_colors, lwd=5, cex=0.9, horiz = TRUE)


if(ressave){
  dev.off()
}




