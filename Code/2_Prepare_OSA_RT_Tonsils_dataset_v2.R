library(dplyr)

calc_association_logistic <- function(pdata, curmodel, covariates, outcomes){
  # curmodel <- mainmodel; pdata <- microbes_species; outcomes <- outcomesFactors
  resP <- as.data.frame(matrix(NA,nrow = length(pdata),ncol = length(outcomes)))
  colnames(resP) <- outcomes; rownames(resP) <- colnames(pdata)
  resB <- as.data.frame(matrix(NA,nrow = length(pdata),ncol = length(outcomes)))
  colnames(resB) <- outcomes; rownames(resB) <- colnames(pdata)
  
  for(outcome in outcomes){ # outcome <- outcomes[1]
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
        # print(tcount)
        tcount = tcount+1
        if(!(tcount %% 100)){print(tcount)}
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


#Load Results
print("Loading Data...")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/MainModel_phenotypes_1117people_12pheno_v1.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_strains_1117people_v1.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_1117people_v1.RData")
# microbes0 <- microbes; mainmodel0 <- mainmodel; pathways0 <- pathways

microbes <- read.table("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Tonsils_HardDecision_Freq_ResultsSummary_SPECIES_cutFreq0_v37_forPandas_AgeSex_v4.txt",
                       header = T,sep="\t")
rownames(microbes) <- microbes$ID; cn <- colnames(microbes); 

mainmodel <- microbes[,c("TotReads","Sex","Age","Condition")] 
microbes <- microbes[,cn[6:length(cn)]] 

#convert bacterial naming to SOL format
cn <- sub(".*?_", "", colnames(microbes))
for(c in 1:length(cn)){ # c <- 640
  tstr <- gsub("Unknown_","Unknown.",cn[c],fixed = T)
  tstr <- gsub("_1_","_",tstr,fixed = T); tstr <- gsub("_2_","_",tstr,fixed = T); 
  tstr <- strsplit(tstr,"_",fixed = T)[[1]]
  if(length(tstr)<7){print(c);print(tstr)}
  nn <- paste0("k__",tstr[1],";p__",tstr[2],";c__",tstr[3],";o__",tstr[4],";f__",tstr[5],";g__",tstr[6],
               ";s__",paste0(tstr[7:length(tstr)],collapse = "_"))
  cn[c] <- gsub("Unknown.","Unknown_",nn,fixed = T)
}
colnames(microbes) <- cn
realnames <- colnames(microbes)

samplesBlanks <- rownames(mainmodel[mainmodel$Condition=="BLANK",])
microbesBlank <- microbes[samplesBlanks,]

# Zero the samples for which there are <10 reads - the abundance threshold
#convert abundances to reads
for(r in rownames(microbes)){ # r <- rownames(microbes)[1]
  tr <- mainmodel[r,"TotReads"]
  microbes[r,] <- microbes[r,]*tr
}
microbes[microbes < 10] <- 0

#Remove Blanks and "RTOSA"
microbes <- microbes[!rownames(microbes) %in% samplesBlanks,]
mainmodel <- mainmodel[!rownames(mainmodel) %in% samplesBlanks,]

samplesRTOSA <- rownames(mainmodel[mainmodel$Condition=="RT+OSA",])
microbes <- microbes[!rownames(microbes) %in% samplesRTOSA,]
mainmodel <- mainmodel[!rownames(mainmodel) %in% samplesRTOSA,]

mainmodel$Sex <- as.factor(mainmodel$Sex)
mainmodel$Condition <- as.factor(mainmodel$Condition)

rn <- rownames(microbes);
mainmodel <- mainmodel[rn,]
cn <- colnames(mainmodel); cn[4] <- "RTOSA"; colnames(mainmodel) <- cn

save(mainmodel, file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_MainModel_phenotypes_129people_v2.RData")
save(microbes, file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_Microbiome_strains_129people_v2.RData")

########################## Read Pathways KO
pathwaysKO <- read.csv("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_KO_pred_metagenome_unstrat_descrip.tsv",
                       header = T,sep="\t")

#converto to SOL format
lkup <- read.csv("/data/linkage/HCHS_SOL/Omics/Microbiome/MS1107 Gut microbiota composition/REP82_ko-enzyme-annotations.txt",
                   header = F,sep="\t")
lkup <- lkup %>% 
  group_by(V1) %>% 
  mutate(V5 = paste0(V5, collapse = " / "))
lkup = as.data.frame(lkup[!duplicated(lkup$V1),])
rownames(lkup) <- lkup$V1

description2 <- lkup[pathwaysKO$function.,"V5"]
pathwaysKO$description <- ifelse(is.na(description2),pathwaysKO$description,description2)

#convert format
pathwaysKO$function. <- paste(pathwaysKO$function.,pathwaysKO$description,sep = ";")
pathwaysKO <- pathwaysKO[ , -which(colnames(pathwaysKO) %in% c("description"))]
pathwaysKO <- as.data.frame(t(pathwaysKO)); colnames(pathwaysKO) <- pathwaysKO[1,]; 
rownames(pathwaysKO) <- gsub(".","/",rownames(pathwaysKO),fixed = T)
pathwaysKO <- pathwaysKO[rn,]
save(pathwaysKO, file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_Microbiome_KO_129people_v2.RData")

########################## Read Pathways KEGG models
pathwaysKEGGModule <- read.csv("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_KEGG_module_abun_unstrat_descrip.tsv",
                       header = T,sep="\t")

#convert format
pathwaysKEGGModule$pathway <- paste(pathwaysKEGGModule$pathway,pathwaysKEGGModule$description,sep = ";")
pathwaysKEGGModule <- pathwaysKEGGModule[ , -which(colnames(pathwaysKEGGModule) %in% c("description"))]
pathwaysKEGGModule <- as.data.frame(t(pathwaysKEGGModule)); colnames(pathwaysKEGGModule) <- pathwaysKEGGModule[1,]; 
rownames(pathwaysKEGGModule) <- gsub(".","/",rownames(pathwaysKEGGModule),fixed = T)
pathwaysKEGGModule <- pathwaysKEGGModule[rn,]
save(pathwaysKEGGModule, file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_Microbiome_KEGGModule_129people_v2.RData")


##############################################################################################
#                                  Analysis                                                  #
##############################################################################################

# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_MainModel_phenotypes_129people_v2.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_Microbiome_strains_129people_v2.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_Microbiome_KO_129people_v2.RData")
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/RTOSA_Microbiome_KEGGModule_129people_v2.RData")

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


#######################################################################################
####################### main model analysis
#######################################################################################

#define models, covariates
covariates <- c("Age","Sex")
outcomesFactors <- c("RTOSA")

###################### species
RTOSA_MainModel_res <- list()
RTOSA_MainModel_res$SpeciesMicrobesRelativeFactor <- calc_association_logistic(microbes_species,mainmodel,covariates,outcomesFactors)
RTOSA_MainModel_res$SpeciesMicrobesRelativeAll <- list()
RTOSA_MainModel_res$SpeciesMicrobesRelativeAll$pval <- cbind(RTOSA_MainModel_res$SpeciesMicrobesRelativeFactor$pval)
RTOSA_MainModel_res$SpeciesMicrobesRelativeAll$beta <- cbind(RTOSA_MainModel_res$SpeciesMicrobesRelativeFactor$beta)

###################### genus
RTOSA_MainModel_res$GenusMicrobesRelativeFactor <- calc_association_logistic(microbes_genus,mainmodel,covariates,outcomesFactors)
RTOSA_MainModel_res$GenusMicrobesRelativeAll <- list()
RTOSA_MainModel_res$GenusMicrobesRelativeAll$pval <- cbind(RTOSA_MainModel_res$GenusMicrobesRelativeFactor$pval)
RTOSA_MainModel_res$GenusMicrobesRelativeAll$beta <- cbind(RTOSA_MainModel_res$GenusMicrobesRelativeFactor$beta)

###################### family
RTOSA_MainModel_res$FamilyMicrobesRelativeFactor <- calc_association_logistic(microbes_family,mainmodel,covariates,outcomesFactors)
RTOSA_MainModel_res$FamilyMicrobesRelativeAll <- list()
RTOSA_MainModel_res$FamilyMicrobesRelativeAll$pval <- cbind(RTOSA_MainModel_res$FamilyMicrobesRelativeFactor$pval)
RTOSA_MainModel_res$FamilyMicrobesRelativeAll$beta <- cbind(RTOSA_MainModel_res$FamilyMicrobesRelativeFactor$beta)

###################### pathways KO
RTOSA_MainModel_res$MicrobesPathwaysKOFactor <- calc_association_logistic(pathwaysKO,mainmodel,covariates,outcomesFactors)
RTOSA_MainModel_res$MicrobesPathwaysKOAll <- list()
RTOSA_MainModel_res$MicrobesPathwaysKOAll$pval <- cbind(RTOSA_MainModel_res$MicrobesPathwaysKOFactor$pval)
RTOSA_MainModel_res$MicrobesPathwaysKOAll$beta <- cbind(RTOSA_MainModel_res$MicrobesPathwaysKOFactor$beta)

###################### pathways KEGGModule
RTOSA_MainModel_res$MicrobesPathwaysKEGGModuleFactor <- calc_association_logistic(pathwaysKEGGModule,mainmodel,covariates,outcomesFactors)
RTOSA_MainModel_res$MicrobesPathwaysKEGGModuleAll <- list()
RTOSA_MainModel_res$MicrobesPathwaysKEGGModuleAll$pval <- cbind(RTOSA_MainModel_res$MicrobesPathwaysKEGGModuleFactor$pval)
RTOSA_MainModel_res$MicrobesPathwaysKEGGModuleAll$beta <- cbind(RTOSA_MainModel_res$MicrobesPathwaysKEGGModuleFactor$beta)


#combine 
RTOSA_MainModel_res$AllBactPval <- rbind(RTOSA_MainModel_res$SpeciesMicrobesRelativeAll$pval,
                                   RTOSA_MainModel_res$GenusMicrobesRelativeAll$pval,
                                   RTOSA_MainModel_res$FamilyMicrobesRelativeAll$pval)
RTOSA_MainModel_res$AllBactBeta <- rbind(RTOSA_MainModel_res$SpeciesMicrobesRelativeAll$beta,
                                   RTOSA_MainModel_res$GenusMicrobesRelativeAll$beta,
                                   RTOSA_MainModel_res$FamilyMicrobesRelativeAll$beta)
# a <- RTOSA_MainModel_res$AllBactPval

RTOSA_MainModel_res$AllFuncPval <- rbind(RTOSA_MainModel_res$MicrobesPathwaysKOAll$pval,
                                         RTOSA_MainModel_res$MicrobesPathwaysKEGGModuleAll$pval)
RTOSA_MainModel_res$AllFuncBeta <- rbind(RTOSA_MainModel_res$MicrobesPathwaysKOAll$beta,
                                            RTOSA_MainModel_res$MicrobesPathwaysKEGGModuleAll$beta)

save(RTOSA_MainModel_res,file = paste0("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/RTOSA_MainModel_v3.RData"))



