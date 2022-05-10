library(R.utils)

addSuffix <- function(model,suffix){
  #model <- MainModel_res$AllBactPval
  cn <- colnames(model); cn <- paste0(cn,suffix)
}


##### Bact Pval #######

AllModels <- list()

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_train4_noranknorm_permutepval_v2.RData")
tmodel <- MainModel_res$AllBactPval; 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_Main")
stopifnot(all(rownames(AllModels)==rownames(tmodel)))
AllModels <- tmodel

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_train4_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModel_res$AllBactPval; 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLife")
stopifnot(all(rownames(AllModels)==rownames(tmodel)))
x <- order(c(1:ncol(AllModels), 1:ncol(tmodel)))
AllModels <- cbind(AllModels,tmodel)[,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModelDiabetes_Lifestyle_Comorbidities_train4_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllBactPval; 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLifeMorb")
stopifnot(all(rownames(AllModels)==rownames(tmodel)))
x <- 1:length(AllModels)
x <- insert(x,ats=seq(3,2*length(tmodel)+2,by=2),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- cbind(AllModels,tmodel)[,x]

#calc min
AllModels$min <- apply(AllModels, 1, FUN = min, na.rm = TRUE)
AllModels <- AllModels[AllModels$min<0.05,]
goodrows <- rownames(AllModels)

write.table(AllModels,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Bacteria_train4_pval_v4.txt",
            sep = "\t", quote = F)

##########################
##### Pathway Pval #######
##########################

AllModels <- list()

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/MainModel_train4_noranknorm_permutepval_v2.RData")
tmodel <- MainModel_res$AllFuncPval; 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_Main")
stopifnot(all(rownames(AllModels)==rownames(tmodel)))
AllModels <- tmodel

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModel_Lifestyle_train4_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModel_res$AllFuncPval; 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLife")
stopifnot(all(rownames(AllModels)==rownames(tmodel)))
x <- order(c(1:ncol(AllModels), 1:ncol(tmodel)))
AllModels <- cbind(AllModels,tmodel)[,x]

load("/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/ExtendedModelDiabetes_Lifestyle_Comorbidities_train4_noranknorm_permutepval_v3.RData")
tmodel <- ExtendedModelDiabetesObesityDiet_res$AllFuncPval; 
colnames(tmodel) <- addSuffix(tmodel,"_permPval_ExtLifeMorb")
stopifnot(all(rownames(AllModels)==rownames(tmodel)))
x <- 1:length(AllModels)
x <- insert(x,ats=seq(3,2*length(tmodel)+2,by=2),values=(length(AllModels)+1):(length(AllModels)+length(tmodel)))
AllModels <- cbind(AllModels,tmodel)[,x]

#calc min
AllModels$min <- apply(AllModels, 1, FUN = min, na.rm = TRUE)
AllModels <- AllModels[AllModels$min<0.05,]
goodrows <- rownames(AllModels)

write.table(AllModels,file = "/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Perm_LifeComorb_Function_train4_pval_v4.txt",
            sep = "\t", quote = F)




