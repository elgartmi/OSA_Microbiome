library(sas7bdat)
library(lubridate)
library(table1)

covariates <- read.sas7bdat("/data/linkage/HCHS_SOL/Omics/Microbiome/HCHS MS1107 Data Transfer/INV2/Datasets/gold_part_derv_inv2.sas7bdat")
medications <- read.sas7bdat("/data/linkage/HCHS_SOL/Omics/Microbiome/HCHS MS1107 Data Transfer/INV2/Datasets/gle_inv2.sas7bdat")
dates <- read.sas7bdat("/data/linkage/HCHS_SOL/Omics/Microbiome/HCHS MS1107 Data Transfer/INV2/Datasets/gol_inv2.sas7bdat")
consent <- read.sas7bdat("/data/linkage/HCHS_SOL/Omics/Microbiome/HCHS MS1107 Data Transfer/INV2/Datasets/gic_inv2.sas7bdat")

clinical <- read.sas7bdat("/data/linkage/HCHS_SOL/Omics/Microbiome/HCHS MS1107 Data Transfer/part_derv_v2_inv3.sas7bdat")
# bb <- read.sas7bdat("/data/linkage/HCHS_SOL/Omics/Microbiome/HCHS MS1107 Data Transfer/lab_derv_v2_inv3.sas7bdat")

microbes <- read.table(file="/data/linkage/HCHS_SOL/Omics/Microbiome/MS1107 Gut microbiota composition/reseq_GOLD_species_absolute_abundance.txt",
                       sep="\t", row.names = 1)
cn <- microbes[1,]; cn <- gsub("R","",cn,fixed = T); cn <- gsub("L","",cn,fixed = T);
colnames(microbes)<-cn; microbes <- microbes[2:dim(microbes)[1],];microbes <- as.data.frame(t(microbes))

covariates$IDn <- as.numeric(covariates$ID)
rownames(covariates) <- covariates$IDn
covariates <- covariates[covariates$BURK_ID!="",]

pathways <- read.table(file="/data/linkage/HCHS_SOL/Omics/Microbiome/MS1107 Gut microbiota composition/reseq_GOLD_KEGG_module_annotated.txt",
                       sep="\t")
cn <- pathways[1,]; 
colnames(pathways)<-cn; pathways <- pathways[2:dim(pathways)[1],]
rownames(pathways) <- paste(pathways$KEGGID,pathways$Level1,pathways$Level2,pathways$Level3,pathways$ModuleID,sep = ";")
pathways <- pathways[,colnames(pathways)[6:length(pathways)]]
cn <- colnames(pathways); cn <- gsub("R","",cn,fixed = T); cn <- gsub("L","",cn,fixed = T);colnames(pathways)<-cn;
pathways <- as.data.frame(t(pathways))

# kegg <- read.table(file="/data/linkage/HCHS_SOL/Omics/Microbiome/MS1107 Gut microbiota composition/reseq_GOLD_KEGG_module_annotated.txt",
#                         sep="\t")
# 

#filter people without consent
consent <- consent[consent$ID %in% covariates$ID,]
consent <- consent[consent$GIC1==1,]
consent <- consent[consent$GIC2==1,]
consent <- consent[consent$GIC3==1,]
consent <- consent[consent$GIC4==1,]
consent <- consent[consent$GIC5==1,]

#filter people that have been taking medications in last 4 months

#antibiotics
medications <- medications[medications$ID %in% consent$ID,]
medications$GLE0A <- as.Date(medications$GLE0A,'1960-01-01')
medications$GLE2A <- as.Date(ISOdate(medications$GLE2AY,medications$GLE2AM,ifelse(is.nan(medications$GLE2AD),1,medications$GLE2AD)))
medications$Interval_Antibiotic <- interval(medications$GLE2A,medications$GLE0A) %/% months(1)

#indigestion
medications$GLE4A <- as.Date(ISOdate(medications$GLE4AY,medications$GLE4AM,ifelse(is.nan(medications$GLE4AD),1,medications$GLE4AD)))
medications$Interval_Indigest <- interval(medications$GLE4A,medications$GLE0A) %/% months(1)

usedata <- medications[((medications$GLE2==2) | (medications$Interval_Antibiotic>3)),]
usedata <- usedata[((usedata$GLE4==2) | (usedata$Interval_Indigest>3)),]

#make file with covariates (age, sex, center, weoght_norm,education, income level, 
#                           smoking, physical activity, alcohol consumption, hypertension,
#                           obesity and diabetes measures ) 
#add outcomes (AHI, AHI_categorical, MinSpO2, AvgSpO2, SpO290)

clinical <- clinical[clinical$ID %in% usedata$ID,]
clinical$IDn <- as.numeric(clinical$ID)
rownames(clinical) <- clinical$IDn

#load OSA outcomes and visit1 data
load("/data/linkage/HCHS_SOL/Projects/2019_genetic_correlations/Phenotypes/phen_raw_20200805.RData")
phenotypes <- phen_raw
rownames(phenotypes) <- phenotypes$HCHS_ID
iids <- which(clinical$IDn %in% phenotypes$HCHS_ID)
clinical <- clinical[iids,]
phenotypes <- phenotypes[as.character(clinical$IDn),]

phenotypes$Weight_V1 <- phenotypes$BMI * (phenotypes$Height/100)^2
clinical$Weight_V2 <- clinical$BMI_V2 * (phenotypes$Height/100)^2

phenotypes2 <- cbind(phenotypes,clinical)
#exclude people with weight change >10% from visit 1 to visit 2
phenotypes2$IncWeight <- (phenotypes2$Weight_V2-phenotypes2$Weight_V1)/phenotypes2$Weight_V2
phenotypes2 <- phenotypes2[phenotypes2$IncWeight < 0.1,]

#add microbial data (species, kegg, enzyme)
# idcv <- read.csv(file="/data/linkage/HCHS_SOL/Omics/Microbiome/MS1107 Gut microbiota composition/ID match.csv")
# rownames(idcv) <- idcv$ID
# idcv <- idcv[phenotypes2$HCHS_ID,]

phenotypes2 <- phenotypes2[phenotypes2$HCHS_ID %in% covariates$IDn,]
covariates <- covariates[as.character(phenotypes2$HCHS_ID),]
phenotypes2$BURK_ID <- covariates$BURK_ID
rownames(phenotypes2) <- phenotypes2$BURK_ID

iids <- which((!is.na(phenotypes2$AHI)) |
                (!is.na(phenotypes2$MinSpO2)) |
                (!is.na(phenotypes2$AvgSpO2)) |
                (!is.na(phenotypes2$SpO290)))
phenotypes2 <- phenotypes2[iids,]
phenotypes2$AHI_C4 <- ifelse(phenotypes2$AHI<5,"no_OSA",
                             ifelse(phenotypes2$AHI<15,"mild_OSA",
                             ifelse(phenotypes2$AHI<30,"moderate_OSA","severe_OSA")))
clinical <- clinical[as.character(phenotypes2$HCHS_ID),]


# after internal review - add diet data
dietdata <- read.csv("/data/linkage/HCHS_SOL/Projects/2021_Puerto_Rican_metabolites/Data/SOL_variables/soldiet_covariates_20210329.csv")
iids <- match(phenotypes2$HCHS_ID,dietdata$ID)
dietdata <- dietdata[iids,]
phenotypes2$mediter_diet_score <- dietdata$mediter_diet_score

allpheno <- cbind(phenotypes2,clinical)


#age, sex, center, weight_norm,education, income level, 
#                           smoking, alcohol consumption, hypertension,
#                           obesity and diabetes measures ) 
#AHI, AHI_categorical, MinSpO2, AvgSpO2, SpO290)
shortphenoN1 <- c("HCHS_ID","scanID","Sex","AGE","BMI","AHI","AHI_C4","MinSpO2","AvgSpO2","SpO290","Height","mediter_diet_score")
shortphenoN2 <- c("CENTER","WEIGHT_NORM_OVERALL_V2",
                  "EDUCATION_C3_V2","INCOME_C3_V2","INCOME_C5_V2",
                  "CURRENT_SMOKER_V2", # "CIGARETTE_PACK_YEARS_C3_V2","CIGARETTE_PACK_YEARS_V2","CIGARETTE_USE_V2","CIGARETTES_YEAR_V2",
                  "ALCOHOL_INTAKE_V2","BMI_V2","BMIGRP_C4_V2","BMIGRP_C6_V2","WAIST_HIP_V2",
                  "HYPERTENSION2_V2","HYPERTENSION_C4_V2","HYPERTENSION2_AHA_V2","HYPERTENSION_AHA_C5_V2","HYPERT_AWARENESS_V2","HYPERT_TREATMENT_V2","HYPERT_CONTROL_V2",
                  "DIABETES_LAB_V2","DIABETES_SELF_V2",
                  "DIABETES3_C4_V2","DIABETES3_INDICATOR_V2","DIABETES3_V2","DIABETES4_C4_V2",
                  "DIABETES4_INDICATOR_V2","DIABETES4_V2","DIABETES5_C4_V2",
                  "DIABETES5_INDICATOR_V2","DIABETES5_V2","DM_TRT_V2","DM3_AWARE_V2",
                  "DM3_CONTROL_V2","DM4_AWARE_V2","DM4_CONTROL_V2",
                  "DM5_AWARE_V2","DM5_CONTROL_V2")

extendedmodel <- cbind(phenotypes2[,shortphenoN1],clinical[,shortphenoN2])

for (c in colnames(extendedmodel)){ #c="Sex"
  td <- as.factor(extendedmodel[,c])
  if(length(levels(td))<8){
    extendedmodel[,c] <- td
  }
}


mainmodelN <- c("HCHS_ID","scanID","Sex","AGE","BMI","AHI","AHI_C4","MinSpO2","AvgSpO2","SpO290",
                "CENTER","WEIGHT_NORM_OVERALL_V2")
mainmodel <- extendedmodel[,mainmodelN] 


rn <- rownames(phenotypes2)[which(rownames(phenotypes2) %in% rownames(microbes))]

allpheno <- allpheno[rn,]
extendedmodel <- extendedmodel[rn,]
mainmodel <- mainmodel[rn,]

#####table1
phenotypes1 <- extendedmodel[,c("AHI","MinSpO2","AvgSpO2","SpO290",
                                "Sex","AGE","CENTER",
                                "EDUCATION_C3_V2", "INCOME_C5_V2", "CURRENT_SMOKER_V2","ALCOHOL_INTAKE_V2", "mediter_diet_score",
                                "HYPERTENSION_C4_V2", "DIABETES4_V2","BMI_V2" ,"WAIST_HIP_V2")]
colnames(phenotypes1) <- c("AHI","MinSpO2","AvgSpO2","SpO290",
                           "Sex","Age","Center",
                           "Education", "Income", "Smoker", "Alcohol","MDS",
                           "Hypertension","Diabetes","BMI","WaistHip")
levels(phenotypes1$Education) <- c("High School", "Below High School", "Above High School","Unknown")
levels(phenotypes1$Income) <- c("< $10,000","$10,001-$20,000","$20,001-$40,000","$40,001-$75,000","> $75,000","Unknown")
levels(phenotypes1$Smoker) <- c("No","Yes")
levels(phenotypes1$Hypertension) <- c("No hypertension","Pre-hypertension","Treated hypertension","Untreated hypertension")
levels(phenotypes1$Center) <- c("Bronx","Chicago","Miami","San Diego")
levels(phenotypes1$Diabetes) <- c("Normal glucose regulation ","Impaired glucose tolerance","Diabetes")

label(phenotypes1$Alcohol) <- "Alcohol Consumption"
label(phenotypes1$Smoker) <- "Current Smoker"
phenotypes1$AHI_C2_at5 <- ifelse(phenotypes1$AHI<5,"no_OSA","OSA")
phenotypes1$AHI_C2_at15 <- ifelse(phenotypes1$AHI<15,"no_OSA","OSA")
phenotypes1$AHI_C2_at30 <- ifelse(phenotypes1$AHI<30,"no_OSA","OSA")
table1(~ AHI + MinSpO2 + AvgSpO2 + SpO290 + Age + Center + Education + Income + Smoker + Alcohol + MDS + Hypertension + Diabetes + BMI + WaistHip + AHI_C2_at5 + AHI_C2_at15 + AHI_C2_at30 | Sex, data=phenotypes1)

###########

save(allpheno,file="/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/all_phenotypes_1117people_436pheno_v2.RData")
save(extendedmodel,file="/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/ExtendedModel_phenotypes_1117people_48pheno_v2.RData")
save(mainmodel,file="/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/MainModel_phenotypes_1117people_12pheno_v2.RData")

microbes <- microbes[rn,]
#remove empty microbes
for(b in colnames(microbes)){
  microbes[,b] <- as.numeric(microbes[,b])
}
cs <- colSums(microbes,na.rm = T)
microbes <- microbes[,which(cs>0)]

save(microbes,file="/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_strains_1117people_v1.RData")

pathways <- pathways[rn,]
#remove empty microbes
for(b in colnames(pathways)){
  pathways[,b] <- as.numeric(pathways[,b])
}
cs <- colSums(pathways,na.rm = T)
pathways <- pathways[,which(cs>0)]

save(pathways,file="/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_1117people_v1.RData")

###################### pathways KO
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_1117people_v1.RData")
# rn <- rownames(pathways)
pathwaysKO <- read.table("/data/linkage/HCHS_SOL/Omics/Microbiome/MS1107 Gut microbiota composition/reseq_GOLD_KEGG.txt",
                         header = T,sep="\t")

lkup <- read.table("/data/linkage/HCHS_SOL/Omics/Microbiome/MS1107 Gut microbiota composition/REP82_ko-enzyme-annotations.txt",
                   header = F,sep="\t")
lkup <- lkup %>% 
  group_by(V1) %>% 
  mutate(V5 = paste0(V5, collapse = " / "))
lkup = as.data.frame(lkup[!duplicated(lkup$V1),])
rownames(lkup) <- lkup$V1

lkup2 <- read.table("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/picrust_out/ko_info.tsv",
                   header = F,sep="\t")
lkup2 <- lkup2 %>% 
  group_by(V1) %>% 
  mutate(V2 = paste0(V2, collapse = " / "))
lkup2 = as.data.frame(lkup2[!duplicated(lkup2$V1),])
rownames(lkup2) <- lkup2$V1

description1 <- lkup[pathwaysKO$X.KEGG.ID,"V5"]
description2 <- lkup2[pathwaysKO$X.KEGG.ID,"V2"]
description <- ifelse(is.na(description1),description2,description1)
description <- ifelse(is.na(description),pathwaysKO$X.KEGG.ID,description)

rownames(pathwaysKO) <- paste0(pathwaysKO$X.KEGG.ID,";",description)
pathwaysKO <- pathwaysKO[,colnames(pathwaysKO)[2:length(pathwaysKO)]]
cn <- colnames(pathwaysKO); cn <- gsub("R","",cn,fixed = T); cn <- gsub("L","",cn,fixed = T);colnames(pathwaysKO)<-cn;
pathwaysKO <- as.data.frame(t(pathwaysKO))
pathwaysKO <- pathwaysKO[rn,]

#remove empty microbes
for(b in colnames(pathwaysKO)){
  pathwaysKO[,b] <- as.numeric(pathwaysKO[,b])
}
cs <- colSums(pathwaysKO,na.rm = T)
pathwaysKO <- pathwaysKO[,which(cs>0)]

save(pathwaysKO,file="/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_PathwaysKO_1117people_v1.RData")

###################### pathways KEGGModule
# load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_Pathways_1117people_v1.RData")
# rn <- rownames(pathways)
pathwaysKEGGModule <- read.table("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/SOL_Kaplan_KEGG_module_abun_unstrat_descrip.tsv",
                         header = T,sep="\t")
#convert format
pathwaysKEGGModule$pathway <- paste(pathwaysKEGGModule$pathway,pathwaysKEGGModule$description,sep = ";")
pathwaysKEGGModule <- pathwaysKEGGModule[ , -which(colnames(pathwaysKEGGModule) %in% c("description"))]
cn <- colnames(pathwaysKEGGModule); cn <- gsub("R","",cn,fixed = T); cn <- gsub("L","",cn,fixed = T);colnames(pathwaysKEGGModule)<-cn;

pathwaysKEGGModule <- as.data.frame(t(pathwaysKEGGModule)); colnames(pathwaysKEGGModule) <- pathwaysKEGGModule[1,]; 
pathwaysKEGGModule <- pathwaysKEGGModule[rn,]

#remove empty microbes
for(b in colnames(pathwaysKEGGModule)){
  pathwaysKEGGModule[,b] <- as.numeric(pathwaysKEGGModule[,b])
}
cs <- colSums(pathwaysKEGGModule,na.rm = T)
pathwaysKEGGModule <- pathwaysKEGGModule[,which(cs>0)]

save(pathwaysKEGGModule,file="/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Data/Microbiome_pathwaysKEGGModule_1117people_v1.RData")




