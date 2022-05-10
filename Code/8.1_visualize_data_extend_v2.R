library(dplyr)
library(eulerr)
library(dendextend)
library(qgraph)
library("psych")
library(ggplot2)

############################################load only significant bacteria from association analysis

#read results GOLD file
GOLDres <- read.csv(file="/data/linkage/HCHS_SOL/Projects/2021_OSA_microbiome/Results/AllModels_Bacteria_pval_v3.txt",sep="\t")
GOLDres <- GOLDres[,!grepl("rnknrm",colnames(GOLDres),fixed = T)]
GOLDres <- GOLDres[,!grepl("_Main",colnames(GOLDres),fixed = T)]

### bypvalue all models
thres <- 0.05

#filter GOLD by bacteria with at least any of the OSA measurements significant (thres) which are also in RTOSA
GOLDres_rn <- c()
GOLDres_rnList <- list()

c0 <- colnames(GOLDres)[grepl("AHI_permPval",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[1]] <- unique(c(rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, unique(c(rownames(a)))))

c0 <- colnames(GOLDres)[grepl("AHI_C2_at30_perm",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[1]] <- unique(c(GOLDres_rnList[[1]],rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, unique(c(rownames(a)))))
GOLDres$AHIallmodels <- rep(FALSE,dim(GOLDres)[1]);GOLDres[GOLDres_rnList[[1]],"AHIallmodels"] <- TRUE

c0 <- colnames(GOLDres)[grepl("MinSpO2_permPval",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[2]] <- unique(c(rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))
GOLDres$MinSpO2allmodels <- rep(FALSE,dim(GOLDres)[1]);GOLDres[GOLDres_rnList[[2]],"MinSpO2allmodels"] <- TRUE

c0 <- colnames(GOLDres)[grepl("AvgSpO2_permPval",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[3]] <- unique(c(rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))
GOLDres$AvgSpO2allmodels <- rep(FALSE,dim(GOLDres)[1]);GOLDres[GOLDres_rnList[[3]],"AvgSpO2allmodels"] <- TRUE

c0 <- colnames(GOLDres)[grepl("SpO290_permPval",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 3, ]; a$Max <- apply(a, 1, max,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[4]] <- unique(c(rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))
GOLDres$SpO290allmodels <- rep(FALSE,dim(GOLDres)[1]);GOLDres[GOLDres_rnList[[4]],"SpO290allmodels"] <- TRUE

c0 <- colnames(GOLDres)[grepl("RTOSA",colnames(GOLDres),fixed = T)]
a <- GOLDres[,c0]; a <- a[rowSums(is.na(a)) < 2, ]; a$Max <- apply(a, 1, min,na.rm = T);a <- a[a$Max < thres,]
GOLDres_rnList[[5]] <- unique(c(rownames(a)))
GOLDres_rn <- unique(c(GOLDres_rn, rownames(a)))
GOLDres$RTOSAallmodels <- rep(FALSE,dim(GOLDres)[1]);GOLDres[GOLDres_rnList[[5]],"RTOSAallmodels"] <- TRUE


GOLDres <- GOLDres[GOLDres_rn,]
outcomespval <- c("RTOSAallmodels","AHIallmodels","MinSpO2allmodels","AvgSpO2allmodels","SpO290allmodels")

# Euler(Venn) plot


GOLDrespval <- GOLDres
GOLDrespval <- GOLDrespval[,outcomespval]

#remove empty
GOLDrespval <- GOLDrespval[rowSums(GOLDrespval)>0,]

#Now need to aggregate for display purposes on family
rn <- rownames(GOLDrespval)
rn2 <- rownames(GOLDrespval)
for(i in 1:length(rn)){ #i=25
  r=rn[i]
  if(grepl("k__",r,fixed = T)){
    tstr <- strsplit(r,";",fixed = T)
    fstr <- paste0(strsplit(tstr[[1]][5],"__")[[1]][2])
    iidx <- which(grepl(fstr,rn,fixed = t))
    if(length(iidx)>0){
      rn[i] <- strsplit(tstr[[1]][5],"__")[[1]][2]
    }else{
      gstr <- paste0("Genus__",strsplit(tstr[[1]][6],"__")[[1]][2])
      iidx <- which(grepl(gstr,rn,fixed = t))
      if(length(iidx)>0){
        rn[i] <- fstr
      }else{
        rn[i] <- strsplit(tstr[[1]][6],"__")[[1]][2]
      }
    }
  }else{
    rn[i] <- strsplit(r,"__")[[1]][2]
  }
}
GOLDrespval$simplename <- rn

GOLDres2 <- GOLDrespval %>% 
  group_by(simplename) %>% 
  summarise(across(everything(), ~max(.x,na.rm = TRUE)))

GOLDres2 <- as.data.frame(GOLDres2)
rn <- GOLDres2$simplename
GOLDres2 <- sapply(GOLDres2[,-1],as.logical)
rownames(GOLDres2) <- rn

fit2 <- euler(GOLDres2, shape = "ellipse")
plot(fit2,lty = 1:3, quantities = T,
     labels = list(font = 4))

fit3 <- euler(GOLDres2[,-1], shape = "ellipse")
plot(fit3,lty = 1:3, quantities = T,
     labels = list(font = 4))



  

