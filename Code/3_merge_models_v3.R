library(R.utils)

addSuffix <- function(model,suffix){
  #model <- MainModel_res$AllBactPval
  cn <- colnames(model); cn <- paste0(cn,suffix)
}


##### Bact Pval #######

AllModels <- list()
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_ANCOMBCANCOM_v3.RData")
tmodel <- MainModel_res$AllBactPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                       "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                       "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                       "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_Main")
AllModels <- tmodel

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_noranknorm_permutepval_v2.RData")
tmodel <- MainModel_res$AllBactPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                       "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                       "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                       "AHI_C2_at30_FDR")] 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_Main")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- order(c(1:ncol(AllModels), 1:ncol(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_ANCOMBCANCOM_v3.RData")
tmodel <- ExtendedModel_res$AllBactPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                           "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                           "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                           "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_ExtLife")
x <- 1:length(AllModels)
x <- insert(x,ats=seq(3,2*length(tmodel)+2,by=2),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModel_res$AllBactPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                           "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                           "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                           "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLife")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- 1:length(AllModels)
x <- insert(x,ats=seq(4,3*length(tmodel)+3,by=3),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_Comorbidities_ANCOMBCANCOM_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllBactPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                                              "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                                              "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                                              "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_ExtLifeMorb")
x <- 1:length(AllModels)
x <- insert(x,ats=seq(5,4*length(tmodel)+4,by=4),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_Comorbidities_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllBactPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                                   "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                                   "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                                   "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLifeMorb")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- 1:length(AllModels)
x <- insert(x,ats=seq(6,5*length(tmodel)+5,by=5),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


#calc min
AllModels$min <- apply(AllModels, 1, FUN = min, na.rm = TRUE)
AllModels <- AllModels[AllModels$min<0.01,]
goodrows <- rownames(AllModels)

write.table(AllModels,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_ANCOM_Perm_LifeComorb_Bacteria_pval_v4.txt",
            sep = "\t", quote = F)


##### Bact beta #######

AllModels <- list()
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_ANCOMBCANCOM_v3.RData")
tmodel <- MainModel_res$AllBactBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                       "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                       "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                       "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_Main")
AllModels <- tmodel

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_noranknorm_permutepval_v2.RData")
tmodel <- MainModel_res$AllBactBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                       "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                       "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                       "AHI_C2_at30_FDR")] 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_Main")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- order(c(1:ncol(AllModels), 1:ncol(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_ANCOMBCANCOM_v3.RData")
tmodel <- ExtendedModel_res$AllBactBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                           "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                           "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                           "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_ExtLife")
x <- 1:length(AllModels)
x <- insert(x,ats=seq(3,2*length(tmodel)+2,by=2),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModel_res$AllBactBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                           "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                           "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                           "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLife")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- 1:length(AllModels)
x <- insert(x,ats=seq(4,3*length(tmodel)+3,by=3),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_Comorbidities_ANCOMBCANCOM_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllBactBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                                              "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                                              "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                                              "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_ExtLifeMorb")
x <- 1:length(AllModels)
x <- insert(x,ats=seq(5,4*length(tmodel)+4,by=4),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_Comorbidities_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllBactBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                                              "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                                              "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                                              "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLifeMorb")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- 1:length(AllModels)
x <- insert(x,ats=seq(6,5*length(tmodel)+5,by=5),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


#calc min
AllModels$min <- apply(AllModels, 1, FUN = min, na.rm = TRUE)
AllModels <- AllModels[AllModels$min<0.01,]
goodrows <- rownames(AllModels)

write.table(AllModels,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_ANCOM_Perm_LifeComorb_Bacteria_beta_v4.txt",
            sep = "\t", quote = F)


##########################
##### Pathway Pval #######
##########################

AllModels <- list()
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_ANCOMBCANCOM_v3.RData")
tmodel <- MainModel_res$AllFuncPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                       "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                       "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                       "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_Main")
AllModels <- tmodel

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_noranknorm_permutepval_v2.RData")
tmodel <- MainModel_res$AllFuncPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                       "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                       "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                       "AHI_C2_at30_FDR")] 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_Main")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- order(c(1:ncol(AllModels), 1:ncol(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_ANCOMBCANCOM_v3.RData")
tmodel <- ExtendedModel_res$AllFuncPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                           "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                           "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                           "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_ExtLife")
x <- 1:length(AllModels)
x <- insert(x,ats=seq(3,2*length(tmodel)+2,by=2),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModel_res$AllFuncPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                           "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                           "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                           "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLife")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- 1:length(AllModels)
x <- insert(x,ats=seq(4,3*length(tmodel)+3,by=3),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_Comorbidities_ANCOMBCANCOM_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllFuncPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                                              "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                                              "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                                              "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_ExtLifeMorb")
x <- 1:length(AllModels)
x <- insert(x,ats=seq(5,4*length(tmodel)+4,by=4),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_Comorbidities_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllFuncPval[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                                              "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                                              "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                                              "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLifeMorb")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- 1:length(AllModels)
x <- insert(x,ats=seq(6,5*length(tmodel)+5,by=5),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


#calc min
AllModels$min <- apply(AllModels, 1, FUN = min, na.rm = TRUE)
AllModels <- AllModels[AllModels$min<0.01,]
goodrows <- rownames(AllModels)

write.table(AllModels,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_ANCOM_Perm_LifeComorb_Function_pval_v4.txt",
            sep = "\t", quote = F)


##### Path beta #######

AllModels <- list()
load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_ANCOMBCANCOM_v3.RData")
tmodel <- MainModel_res$AllFuncBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                       "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                       "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                       "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_Main")
AllModels <- tmodel

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_noranknorm_permutepval_v2.RData")
tmodel <- MainModel_res$AllFuncBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                       "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                       "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                       "AHI_C2_at30_FDR")] 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_Main")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- order(c(1:ncol(AllModels), 1:ncol(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_ANCOMBCANCOM_v3.RData")
tmodel <- ExtendedModel_res$AllFuncBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                           "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                           "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                           "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_ExtLife")
x <- 1:length(AllModels)
x <- insert(x,ats=seq(3,2*length(tmodel)+2,by=2),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModel_res$AllFuncBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                           "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                           "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                           "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLife")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- 1:length(AllModels)
x <- insert(x,ats=seq(4,3*length(tmodel)+3,by=3),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_Comorbidities_ANCOMBCANCOM_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllFuncBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                                              "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                                              "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                                              "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_ANCOM_ExtLifeMorb")
x <- 1:length(AllModels)
x <- insert(x,ats=seq(5,4*length(tmodel)+4,by=4),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_Comorbidities_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllFuncBeta[,c("AHI","MinSpO2","AvgSpO2","SpO290","AHI_FDR",
                                                              "MinSpO2_FDR","AvgSpO2_FDR","SpO290_FDR","AHI_C2_at5","AHI_C2_at15",
                                                              "AHI_C2_at30","AHI_C2_at5_FDR","AHI_C2_at15_FDR",
                                                              "AHI_C2_at30_FDR")]
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLifeMorb")
rn <- rownames(tmodel); rn <- gsub("Genus__","",rn); rn <- gsub("Family__","",rn);rownames(tmodel) <- rn
x <- 1:length(AllModels)
x <- insert(x,ats=seq(6,5*length(tmodel)+5,by=5),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- merge(AllModels,tmodel,by="row.names",all = T)
rownames(AllModels) <- AllModels$Row.names
cn <- colnames(AllModels)[2:length(colnames(AllModels))]; 
AllModels <- AllModels[,cn][,x]


#calc min
AllModels$min <- apply(AllModels, 1, FUN = min, na.rm = TRUE)
AllModels <- AllModels[AllModels$min<0.01,]
goodrows <- rownames(AllModels)

write.table(AllModels,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_ANCOM_Perm_LifeComorb_Function_beta_v4.txt",
            sep = "\t", quote = F)


