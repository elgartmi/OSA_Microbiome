library(lmPerm)
library(dplyr)

DORANKNORM <- F
PERMUTEPVAL <- T
# DORANKNORM <- T 
# PERMUTEPVAL <- F

if(DORANKNORM){
  rnstr <- "ranknorm"  
}else{
  rnstr <- "noranknorm"  
}

if(PERMUTEPVAL){
  rnstr <- paste0(rnstr,"_permutepval")
}else{
  rnstr <- paste0(rnstr,"_NOpermutepval")
}


calc_association_linear <- function(pdata, curmodel, covariates, outcomes, ranknorm = T, permutepval = F){
  # curmodel <- mainmodel; pdata <- microbes_species; permutepval = PERMUTEPVAL; ranknorm = DORANKNORM
  # curmodel <- mainmodel; pdata <- microbes_genus
  # curmodel <- extendedmodel; pdata <- microbes_species
  # curmodel <- mainmodel; pdata <- microbes_genus
  resP <- as.data.frame(matrix(NA,nrow = length(pdata),ncol = length(outcomes)))
  colnames(resP) <- outcomes; rownames(resP) <- colnames(pdata)
  resB <- as.data.frame(matrix(NA,nrow = length(pdata),ncol = length(outcomes)))
  colnames(resB) <- outcomes; rownames(resB) <- colnames(pdata)
  
  for(outcome in outcomes){ # outcome <- outcomes[1]
    print(outcome)
    #ranknorm
    x <-curmodel[,outcome]
    if(ranknorm==T){
      x <- scale(qnorm(rank(x, ties.method = "random",na.last="keep")/(length(x)+1)))
    }
    tcount=0
    for(b in colnames(pdata)){ # b=colnames(pdata)[1]
      tdata <- curmodel[,c(outcome,covariates)]
      tdata[,1] <- x
      y <- as.numeric(pdata[,b])
      tdata[,'bact'] <- y
      model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(covariates,"bact"),collapse="+")))
      if(permutepval){
        mod <- tryCatch(
          {
            log <- capture.output({
              a<- lmp(model.formula, data = tdata,na.action="na.exclude",
                      center=FALSE,contrasts=list(Sex="contr.treatment",
                      CENTER="contr.treatment"),singular.ok=T) 
            })
            a
          },
          error=function(cond) {
            mod=list()
          }
        )
        pvpos <- 3
      }else{
        mod <- lm(model.formula, data = tdata,na.action="na.exclude") #summary(mod)
        pvpos <- 4
      }
      # r2 <- resid(mod,na.action="na.exclude")
      # shapiro.test(r2)
      if(length(mod)<1){
        # print("Fit failed")
        bactpval <- NA
        bactbeta <- NA
      }else if("bact" %in% rownames(summary(mod)$coefficients)){
        if((permutepval==T)&(summary(mod)$coefficients["bact",pvpos]==0)){
          log <- capture.output({
            mod<- lmp(model.formula, data = tdata,na.action="na.exclude",
                    center=FALSE,contrasts=list(Sex="contr.treatment",
                                                CENTER="contr.treatment"),singular.ok=T,
                    perm="Prob", Ca=0, maxIter=100000) 
          })
        }
        bactpval <- summary(mod)$coefficients["bact",pvpos]
        bactbeta <- summary(mod)$coefficients["bact",1]
      }else{
        bactpval <- 1
        bactbeta <- 0
      }
      resP[b,outcome] <- bactpval
      resB[b,outcome] <- bactbeta
      if(!is.na(bactpval)){
        if(bactpval<0.001){
          print(paste(outcome,b,bactpval))
        }
      }
      tcount = tcount+1
      if(!(tcount %% 100)){print(tcount)}
    }
  }
  
  for(outcome in outcomes){ # outcome <- "AHI"
    p <- resP[,outcome]
    pa <- p.adjust(p,"BH")
    resP[,paste0(outcome,"_FDR")] <- pa
  }
  
  res <- list();res$pval <- resP; res$beta <- resB
  return(res)
}

calc_association_logistic <- function(pdata, curmodel, covariates, outcomes){
  # curmodel <- mainmodel; pdata <- microbes; outcomes <- outcomesFactors
  resP <- as.data.frame(matrix(NA,nrow = length(pdata),ncol = length(outcomes)))
  colnames(resP) <- outcomes; rownames(resP) <- colnames(pdata)
  resB <- as.data.frame(matrix(NA,nrow = length(pdata),ncol = length(outcomes)))
  colnames(resB) <- outcomes; rownames(resB) <- colnames(pdata)
  
  for(outcome in outcomes){ # outcome <- outcomes[3]
    print(outcome)
    tcount=0
    for(b in colnames(pdata)){ # b=colnames(pdata)[1]
      tdata <- curmodel[,c(outcome,covariates)]
      tdata[,'bact'] <- as.numeric(pdata[,b])
      tdata <- tdata[!is.na(tdata[,outcome]),]
      tdata[,outcome] <- as.factor(tdata[,outcome])
      if(sum(tdata[,'bact'],na.rm = T)==0){
        resP[b,outcome] <- 1
        resB[b,outcome] <- 0
        next
      }
      model.formula <- as.formula(paste(paste(outcome," ~ "), paste0(c(covariates,"bact"),collapse="+")))
      
      mod <- glm(model.formula, data = tdata,family=binomial, na.action="na.exclude") #summary(mod)
      # r2 <- resid(mod,na.action="na.exclude")
      # shapiro.test(r2)
      if("bact" %in% rownames(summary(mod)$coefficients)){
        bactpval <- summary(mod)$coefficients["bact",4] 
        bactbeta <- summary(mod)$coefficients["bact",1]
      }else{
        bactpval <- 1
        bactbeta <- 0
      }
      resP[b,outcome] <- bactpval
      resB[b,outcome] <- bactbeta
      if(bactpval<0.001){
        print(paste(outcome,b,bactpval))
      }
      tcount = tcount+1
      if(!(tcount %% 100)){print(tcount)}
    }
  }
  
  for(outcome in outcomes){ # outcome <- "AHI"
    p <- resP[,outcome]
    pa <- p.adjust(p,"BH")
    resP[,paste0(outcome,"_FDR")] <- pa
  }
  
  res <- list();res$pval <- resP; res$beta <- resB
  return(res)
}


#################################### Load Data ##########################################################

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/MainModel_phenotypes_893people_12pheno_train4_v1.RData")
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/ExtendedModel_phenotypes_1117people_48pheno_v2.RData")
extendedmodel_train <- extendedmodel[rownames(mainmodel_train),]
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_strains_893people_train4_v1.RData")
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_893people_train4_v1.RData")
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKO_893people_train4_v1.RData")
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKEGGModule_893people_train4_v1.RData")

mainmodel <- mainmodel_train
extendedmodel <- extendedmodel_train
microbes <- microbes_train
pathways <- pathways_train
pathwaysKO <- pathwaysKO_train
pathwaysKEGGModule <- pathwaysKEGGModule_train



################################## generate new split train/test ############################3

# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/MainModel_phenotypes_1117people_12pheno_v1.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/ExtendedModel_phenotypes_1117people_47pheno_v1.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_strains_1117people_v1.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_1117people_v1.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKO_1117people_v1.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_pathwaysKEGGModule_1117people_v1.RData")
# 
# #option1 - just divide (8 AHI cases)
# iidxtrain <- sample(893)  #nrow(mainmodel)*0.8
# 
# # #option2 - make sure we have 20 AHI cases
# # iidxOSA <- which(mainmodel$AHI>30)
# # iidxnoOSA <- which(mainmodel$AHI<=30)
# # iidxtrain <- c(sample(iidxOSA,23),sample(iidxnoOSA,893-23))
# 
# length(which(mainmodel$AHI>30))
# length(which(mainmodel$AHI[iidxtrain]>30))
# length(which(mainmodel$AHI[-iidxtrain]>30))

# mainmodel_train <- mainmodel[iidxtrain,]
# mainmodel_test <- mainmodel[-iidxtrain,]
# save(mainmodel_train,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/MainModel_phenotypes_893people_12pheno_train4_v1.RData")
# save(mainmodel_test,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/MainModel_phenotypes_224people_12pheno_test4_v1.RData")
# 
# extendedmodel_train <- extendedmodel[iidxtrain,]
# extendedmodel_test <- extendedmodel[-iidxtrain,]
# save(extendedmodel_train,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/ExtendedModel_phenotypes_893people_47pheno_train4_v1.RData")
# save(extendedmodel_test,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/ExtendedModel_phenotypes_224people_47pheno_test4_v1.RData")
# 
# microbes_train <- microbes[iidxtrain,]
# microbes_test <- microbes[-iidxtrain,]
# save(microbes_train,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_strains_893people_train4_v1.RData")
# save(microbes_test,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_strains_224people_test4_v1.RData")
# 
# pathways_train <- pathways[iidxtrain,]
# pathways_test <- pathways[-iidxtrain,]
# save(pathways_train,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_893people_train4_v1.RData")
# save(pathways_test,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_224people_test4_v1.RData")
# 
# pathwaysKO_train <- pathwaysKO[iidxtrain,]
# pathwaysKO_test <- pathwaysKO[-iidxtrain,]
# save(pathwaysKO_train,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKO_893people_train4_v1.RData")
# save(pathwaysKO_test,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKO_224people_test4_v1.RData")
# 
# pathwaysKEGGModule_train <- pathwaysKEGGModule[iidxtrain,]
# pathwaysKEGGModule_test <- pathwaysKEGGModule[-iidxtrain,]
# save(pathwaysKEGGModule_train,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKEGGModule_893people_train4_v1.RData")
# save(pathwaysKEGGModule_test,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKEGGModule_224people_test4_v1.RData")

##############################################################################################
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


outcomes <- c("AHI","MinSpO2","AvgSpO2","SpO290")
outcomesFactors <- c("AHI_C2_at5","AHI_C2_at15","AHI_C2_at30","AHI_C2_at5_15","AHI_C2_at5_30")

#######################################################################################
####################### main model
#######################################################################################
if(FALSE){
  #define models, covariates
  covariates <- c("AGE","Sex","CENTER","WEIGHT_NORM_OVERALL_V2")
  print("Main Model")
  print(covariates)  
  
  #define categories for logistic regression
  mainmodel$AHI_C2_at5 <- ifelse(mainmodel$AHI<5,"no_OSA","OSA")
  mainmodel$AHI_C2_at15 <- ifelse(mainmodel$AHI<15,"no_OSA","OSA")
  mainmodel$AHI_C2_at30 <- ifelse(mainmodel$AHI<30,"no_OSA","OSA")
  mainmodel$AHI_C2_at5_15 <- ifelse(mainmodel$AHI<5,"no_OSA",ifelse(mainmodel$AHI>15,"OSA",NA))
  mainmodel$AHI_C2_at5_30 <- ifelse(mainmodel$AHI<5,"no_OSA",ifelse(mainmodel$AHI>30,"OSA",NA))
  
  ###################### species
  MainModel_res <- list()
  MainModel_res$SpeciesMicrobesRelative <- calc_association_linear(microbes_species,mainmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  MainModel_res$SpeciesMicrobesRelativeFactor <- calc_association_logistic(microbes_species,mainmodel,covariates,outcomesFactors)
  MainModel_res$SpeciesMicrobesRelativeAll <- list()
  MainModel_res$SpeciesMicrobesRelativeAll$pval <- cbind(MainModel_res$SpeciesMicrobesRelative$pval,
                                                    MainModel_res$SpeciesMicrobesRelativeFactor$pval)
  MainModel_res$SpeciesMicrobesRelativeAll$beta <- cbind(MainModel_res$SpeciesMicrobesRelative$beta,
                                                    MainModel_res$SpeciesMicrobesRelativeFactor$beta)
  
  ###################### genus
  MainModel_res$GenusMicrobesRelative <- calc_association_linear(microbes_genus,mainmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  MainModel_res$GenusMicrobesRelativeFactor <- calc_association_logistic(microbes_genus,mainmodel,covariates,outcomesFactors)
  MainModel_res$GenusMicrobesRelativeAll <- list()
  MainModel_res$GenusMicrobesRelativeAll$pval <- cbind(MainModel_res$GenusMicrobesRelative$pval,
                                                       MainModel_res$GenusMicrobesRelativeFactor$pval)
  MainModel_res$GenusMicrobesRelativeAll$beta <- cbind(MainModel_res$GenusMicrobesRelative$beta,
                                                       MainModel_res$GenusMicrobesRelativeFactor$beta)
  
  ###################### family
  MainModel_res$FamilyMicrobesRelative <- calc_association_linear(microbes_family,mainmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  MainModel_res$FamilyMicrobesRelativeFactor <- calc_association_logistic(microbes_family,mainmodel,covariates,outcomesFactors)
  MainModel_res$FamilyMicrobesRelativeAll <- list()
  MainModel_res$FamilyMicrobesRelativeAll$pval <- cbind(MainModel_res$FamilyMicrobesRelative$pval,
                                                        MainModel_res$FamilyMicrobesRelativeFactor$pval)
  MainModel_res$FamilyMicrobesRelativeAll$beta <- cbind(MainModel_res$FamilyMicrobesRelative$beta,
                                                        MainModel_res$FamilyMicrobesRelativeFactor$beta)
  
  ###################### pathways
  MainModel_res$MicrobesPathways <- calc_association_linear(pathways,mainmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  MainModel_res$MicrobesPathwaysFactor <- calc_association_logistic(pathways,mainmodel,covariates,outcomesFactors)
  MainModel_res$MicrobesPathwaysAll <- list()
  MainModel_res$MicrobesPathwaysAll$pval <- cbind(MainModel_res$MicrobesPathways$pval,
                                                  MainModel_res$MicrobesPathwaysFactor$pval)
  MainModel_res$MicrobesPathwaysAll$beta <- cbind(MainModel_res$MicrobesPathways$beta,
                                                  MainModel_res$MicrobesPathwaysFactor$beta)
  
  MainModel_res$MicrobesPathwaysKO <- calc_association_linear(pathwaysKO,mainmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  MainModel_res$MicrobesPathwaysKOFactor <- calc_association_logistic(pathwaysKO,mainmodel,covariates,outcomesFactors)
  MainModel_res$MicrobesPathwaysKOAll <- list()
  MainModel_res$MicrobesPathwaysKOAll$pval <- cbind(MainModel_res$MicrobesPathwaysKO$pval,
                                                    MainModel_res$MicrobesPathwaysKOFactor$pval)
  MainModel_res$MicrobesPathwaysKOAll$beta <- cbind(MainModel_res$MicrobesPathwaysKO$beta,
                                                    MainModel_res$MicrobesPathwaysKOFactor$beta)
  
  MainModel_res$MicrobesPathwaysKEGGModule <- calc_association_linear(pathwaysKEGGModule,mainmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  MainModel_res$MicrobesPathwaysKEGGModuleFactor <- calc_association_logistic(pathwaysKEGGModule,mainmodel,covariates,outcomesFactors)
  MainModel_res$MicrobesPathwaysKEGGModuleAll <- list()
  MainModel_res$MicrobesPathwaysKEGGModuleAll$pval <- cbind(MainModel_res$MicrobesPathwaysKEGGModule$pval,
                                                    MainModel_res$MicrobesPathwaysKEGGModuleFactor$pval)
  MainModel_res$MicrobesPathwaysKEGGModuleAll$beta <- cbind(MainModel_res$MicrobesPathwaysKEGGModule$beta,
                                                    MainModel_res$MicrobesPathwaysKEGGModuleFactor$beta)
  
  
  #combine 
  MainModel_res$AllBactPval <- rbind(MainModel_res$SpeciesMicrobesRelativeAll$pval,
                                     MainModel_res$GenusMicrobesRelativeAll$pval,
                                     MainModel_res$FamilyMicrobesRelativeAll$pval)
  MainModel_res$AllBactBeta <- rbind(MainModel_res$SpeciesMicrobesRelativeAll$beta,
                                     MainModel_res$GenusMicrobesRelativeAll$beta,
                                     MainModel_res$FamilyMicrobesRelativeAll$beta)
  
  MainModel_res$AllFuncPval <- rbind(MainModel_res$MicrobesPathwaysAll$pval,
                                     MainModel_res$MicrobesPathwaysKOAll$pval,
                                     MainModel_res$MicrobesPathwaysKEGGModuleAll$pval)
  MainModel_res$AllFuncBeta <- rbind(MainModel_res$MicrobesPathwaysAll$beta,
                                     MainModel_res$MicrobesPathwaysKOAll$beta,
                                     MainModel_res$MicrobesPathwaysKEGGModuleAll$beta)
  
  
  save(MainModel_res,file = paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_train4_",rnstr,"_v2.RData"))

}

#######################################################################################
####################### Extendedmodel - Lifestyle
#######################################################################################

if(TRUE){
  #define categories for logistic regression
  extendedmodel$AHI_C2_at5 <- ifelse(extendedmodel$AHI<5,"no_OSA","OSA")
  extendedmodel$AHI_C2_at15 <- ifelse(extendedmodel$AHI<15,"no_OSA","OSA")
  extendedmodel$AHI_C2_at30 <- ifelse(extendedmodel$AHI<30,"no_OSA","OSA")
  extendedmodel$AHI_C2_at5_15 <- ifelse(extendedmodel$AHI<5,"no_OSA",ifelse(extendedmodel$AHI>15,"OSA",NA))
  extendedmodel$AHI_C2_at5_30 <- ifelse(extendedmodel$AHI<5,"no_OSA",ifelse(extendedmodel$AHI>30,"OSA",NA))
  
  #define models, covariates
  covariates <- c("AGE","Sex","CENTER","WEIGHT_NORM_OVERALL_V2",
                  "EDUCATION_C3_V2",
                  "INCOME_C5_V2",                               # "INCOME_C3_V2"
                  "CURRENT_SMOKER_V2","ALCOHOL_INTAKE_V2", "mediter_diet_score"
  )
  
  print("Extendedmodel - Lifestyle")
  print(covariates)  
  
  
  ###################### species
  ExtendedModel_res <- list()
  ExtendedModel_res$SpeciesMicrobesRelative <- calc_association_linear(microbes_species,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModel_res$SpeciesMicrobesRelativeFactor <- calc_association_logistic(microbes_species,extendedmodel,covariates,outcomesFactors)
  ExtendedModel_res$SpeciesMicrobesRelativeAll <- list()
  ExtendedModel_res$SpeciesMicrobesRelativeAll$pval <- cbind(ExtendedModel_res$SpeciesMicrobesRelative$pval,
                                                         ExtendedModel_res$SpeciesMicrobesRelativeFactor$pval)
  ExtendedModel_res$SpeciesMicrobesRelativeAll$beta <- cbind(ExtendedModel_res$SpeciesMicrobesRelative$beta,
                                                         ExtendedModel_res$SpeciesMicrobesRelativeFactor$beta)
  
  ###################### genus
  ExtendedModel_res$GenusMicrobesRelative <- calc_association_linear(microbes_genus,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModel_res$GenusMicrobesRelativeFactor <- calc_association_logistic(microbes_genus,extendedmodel,covariates,outcomesFactors)
  ExtendedModel_res$GenusMicrobesRelativeAll <- list()
  ExtendedModel_res$GenusMicrobesRelativeAll$pval <- cbind(ExtendedModel_res$GenusMicrobesRelative$pval,
                                                       ExtendedModel_res$GenusMicrobesRelativeFactor$pval)
  ExtendedModel_res$GenusMicrobesRelativeAll$beta <- cbind(ExtendedModel_res$GenusMicrobesRelative$beta,
                                                       ExtendedModel_res$GenusMicrobesRelativeFactor$beta)
  
  ###################### family
  ExtendedModel_res$FamilyMicrobesRelative <- calc_association_linear(microbes_family,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModel_res$FamilyMicrobesRelativeFactor <- calc_association_logistic(microbes_family,extendedmodel,covariates,outcomesFactors)
  ExtendedModel_res$FamilyMicrobesRelativeAll <- list()
  ExtendedModel_res$FamilyMicrobesRelativeAll$pval <- cbind(ExtendedModel_res$FamilyMicrobesRelative$pval,
                                                        ExtendedModel_res$FamilyMicrobesRelativeFactor$pval)
  ExtendedModel_res$FamilyMicrobesRelativeAll$beta <- cbind(ExtendedModel_res$FamilyMicrobesRelative$beta,
                                                        ExtendedModel_res$FamilyMicrobesRelativeFactor$beta)
  
  ###################### pathways
  ExtendedModel_res$MicrobesPathways <- calc_association_linear(pathways,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModel_res$MicrobesPathwaysFactor <- calc_association_logistic(pathways,extendedmodel,covariates,outcomesFactors)
  ExtendedModel_res$MicrobesPathwaysAll <- list()
  ExtendedModel_res$MicrobesPathwaysAll$pval <- cbind(ExtendedModel_res$MicrobesPathways$pval,
                                                  ExtendedModel_res$MicrobesPathwaysFactor$pval)
  ExtendedModel_res$MicrobesPathwaysAll$beta <- cbind(ExtendedModel_res$MicrobesPathways$beta,
                                                  ExtendedModel_res$MicrobesPathwaysFactor$beta)
  
  ExtendedModel_res$MicrobesPathwaysKO <- calc_association_linear(pathwaysKO,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModel_res$MicrobesPathwaysKOFactor <- calc_association_logistic(pathwaysKO,extendedmodel,covariates,outcomesFactors)
  ExtendedModel_res$MicrobesPathwaysKOAll <- list()
  ExtendedModel_res$MicrobesPathwaysKOAll$pval <- cbind(ExtendedModel_res$MicrobesPathwaysKO$pval,
                                                    ExtendedModel_res$MicrobesPathwaysKOFactor$pval)
  ExtendedModel_res$MicrobesPathwaysKOAll$beta <- cbind(ExtendedModel_res$MicrobesPathwaysKO$beta,
                                                    ExtendedModel_res$MicrobesPathwaysKOFactor$beta)
  
  ExtendedModel_res$MicrobesPathwaysKEGGModule <- calc_association_linear(pathwaysKEGGModule,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModel_res$MicrobesPathwaysKEGGModuleFactor <- calc_association_logistic(pathwaysKEGGModule,extendedmodel,covariates,outcomesFactors)
  ExtendedModel_res$MicrobesPathwaysKEGGModuleAll <- list()
  ExtendedModel_res$MicrobesPathwaysKEGGModuleAll$pval <- cbind(ExtendedModel_res$MicrobesPathwaysKEGGModule$pval,
                                                            ExtendedModel_res$MicrobesPathwaysKEGGModuleFactor$pval)
  ExtendedModel_res$MicrobesPathwaysKEGGModuleAll$beta <- cbind(ExtendedModel_res$MicrobesPathwaysKEGGModule$beta,
                                                            ExtendedModel_res$MicrobesPathwaysKEGGModuleFactor$beta)
  
  #combine 
  ExtendedModel_res$AllBactPval <- rbind(ExtendedModel_res$SpeciesMicrobesRelativeAll$pval,
                                     ExtendedModel_res$GenusMicrobesRelativeAll$pval,
                                     ExtendedModel_res$FamilyMicrobesRelativeAll$pval)
  ExtendedModel_res$AllBactBeta <- rbind(ExtendedModel_res$SpeciesMicrobesRelativeAll$beta,
                                     ExtendedModel_res$GenusMicrobesRelativeAll$beta,
                                     ExtendedModel_res$FamilyMicrobesRelativeAll$beta)
  
  ExtendedModel_res$AllFuncPval <- rbind(ExtendedModel_res$MicrobesPathwaysAll$pval,
                                     ExtendedModel_res$MicrobesPathwaysKOAll$pval,
                                     ExtendedModel_res$MicrobesPathwaysKEGGModuleAll$pval)
  ExtendedModel_res$AllFuncBeta <- rbind(ExtendedModel_res$MicrobesPathwaysAll$beta,
                                     ExtendedModel_res$MicrobesPathwaysKOAll$beta,
                                     ExtendedModel_res$MicrobesPathwaysKEGGModuleAll$beta)
  
  
  save(ExtendedModel_res,file = paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_train4_",rnstr,"_v3.RData"))

}

#######################################################################################
####################### extended model - Lifestyle + Comorbidities (Diabetes, Obesity, Hypertension)
#######################################################################################
if(TRUE){
  #define models, covariates
  covariates <- c("AGE","Sex","CENTER","WEIGHT_NORM_OVERALL_V2",
                  "EDUCATION_C3_V2",
                  "INCOME_C5_V2",                               
                  "CURRENT_SMOKER_V2","ALCOHOL_INTAKE_V2",
                  "HYPERTENSION_C4_V2",
                  "DIABETES4_V2",
                  "BMI_V2","WAIST_HIP_V2","mediter_diet_score"
  )
  print("Extendedmodel - Lifestyle + Comorbidities (Diabetes, Obesity, Hypertension)")
  print(covariates)  
  
  
  ###################### species
  ExtendedModelDiabetesObesityDiet_res <- list()
  ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelative <- calc_association_linear(microbes_species,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelativeFactor <- calc_association_logistic(microbes_species,extendedmodel,covariates,outcomesFactors)
  ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelativeAll <- list()
  ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelativeAll$pval <- cbind(ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelative$pval,
                                                                                ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelativeFactor$pval)
  ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelativeAll$beta <- cbind(ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelative$beta,
                                                                                ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelativeFactor$beta)
  
  ###################### genus
  ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelative <- calc_association_linear(microbes_genus,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelativeFactor <- calc_association_logistic(microbes_genus,extendedmodel,covariates,outcomesFactors)
  ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelativeAll <- list()
  ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelativeAll$pval <- cbind(ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelative$pval,
                                                                              ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelativeFactor$pval)
  ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelativeAll$beta <- cbind(ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelative$beta,
                                                                              ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelativeFactor$beta)
  
  ###################### family
  ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelative <- calc_association_linear(microbes_family,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelativeFactor <- calc_association_logistic(microbes_family,extendedmodel,covariates,outcomesFactors)
  ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelativeAll <- list()
  ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelativeAll$pval <- cbind(ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelative$pval,
                                                                               ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelativeFactor$pval)
  ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelativeAll$beta <- cbind(ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelative$beta,
                                                                               ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelativeFactor$beta)
  
  ###################### pathways
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathways <- calc_association_linear(pathways,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysFactor <- calc_association_logistic(pathways,extendedmodel,covariates,outcomesFactors)
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysAll <- list()
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysAll$pval <- cbind(ExtendedModelDiabetesObesityDiet_res$MicrobesPathways$pval,
                                                                         ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysFactor$pval)
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysAll$beta <- cbind(ExtendedModelDiabetesObesityDiet_res$MicrobesPathways$beta,
                                                                         ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysFactor$beta)
  
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKO <- calc_association_linear(pathwaysKO,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKOFactor <- calc_association_logistic(pathwaysKO,extendedmodel,covariates,outcomesFactors)
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKOAll <- list()
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKOAll$pval <- cbind(ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKO$pval,
                                                                           ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKOFactor$pval)
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKOAll$beta <- cbind(ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKO$beta,
                                                                           ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKOFactor$beta)
  
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModule <- calc_association_linear(pathwaysKEGGModule,extendedmodel,covariates,outcomes,DORANKNORM,PERMUTEPVAL)
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModuleFactor <- calc_association_logistic(pathwaysKEGGModule,extendedmodel,covariates,outcomesFactors)
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModuleAll <- list()
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModuleAll$pval <- cbind(ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModule$pval,
                                                                                   ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModuleFactor$pval)
  ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModuleAll$beta <- cbind(ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModule$beta,
                                                                                   ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModuleFactor$beta)
  
  
  #combine 
  ExtendedModelDiabetesObesityDiet_res$AllBactPval <- rbind(ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelativeAll$pval,
                                                            ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelativeAll$pval,
                                                            ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelativeAll$pval)
  ExtendedModelDiabetesObesityDiet_res$AllBactBeta <- rbind(ExtendedModelDiabetesObesityDiet_res$SpeciesMicrobesRelativeAll$beta,
                                                            ExtendedModelDiabetesObesityDiet_res$GenusMicrobesRelativeAll$beta,
                                                            ExtendedModelDiabetesObesityDiet_res$FamilyMicrobesRelativeAll$beta)
  
  
  ExtendedModelDiabetesObesityDiet_res$AllFuncPval <- rbind(ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysAll$pval,
                                                            ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKOAll$pval,
                                                            ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModuleAll$pval)
  ExtendedModelDiabetesObesityDiet_res$AllFuncBeta <- rbind(ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysAll$beta,
                                                            ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKOAll$beta,
                                                            ExtendedModelDiabetesObesityDiet_res$MicrobesPathwaysKEGGModuleAll$beta)
  
  
  
  save(ExtendedModelDiabetesObesityDiet_res,file = paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModelDiabetes_Lifestyle_Comorbidities_train4_",rnstr,"_v3.RData"))

}

