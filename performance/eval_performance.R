# codes to evaluate performance of the different prediction methods compared

options(scipen = 9999)
#library(reshape)
#library(tidyverse)
#library(emmeans)
#library(dplyr)
#library(reshape2)

# 1. genome-wide performance ####
truecols=c("Recent Base Population","Intermediate Base Population","Distant Base Population")
testcols= c("IBD_Haplo15c", "IBD_Haplo9c","GIBDLD","LocalNgsRelate","TRUFFLE",
          "ZooRoH-125", "ZooRoH-25","ZooRoH-1R","PLINK ROH",
          "phasedibd","hap-IBD","GERMLINE", "Refined IBD",
          "UNI","GRM", "Pedigree")
load('../merge/global_sim_mediumdensity.RData')

# a. Correlations between predicted and reference genome-wide inbreeding
corglob <- cor( data[,testcols], data[,truecols] )

# b. Average genome-wide predicted and reference inbreeding 
ff = data[,4:ncol(data)]
meanF <- colMeans(ff[sapply(ff, is.numeric)], na.rm = TRUE)

# 2. locus-specific performance ####
truecols=c("Recent Base Population","Intermediate Base Population","Distant Base Population")
testcols= c("IBD_Haplo15c", "IBD_Haplo9c","GIBDLD","LocalNgsRelate","TRUFFLE",
          "ZooRoH-125", "ZooRoH-25","ZooRoH-1R","PLINK ROH",
          "phasedibd","hap-IBD","GERMLINE", "Refined IBD")
load('../merge/local_sim_mediumdensity.RData')

# a. Correlations between predicted and reference locus-specific inbreeding 
corloc <- cor( data[,testcols], data[,truecols] )

# for b, c and d sections below, when using real data,
# use the fields corresponding to the estimated HBD status ('Recent Stat', 'Intermediate Stat' and 'Distant Stat')
# they are codes (0 for nonHBD and 1 for HBD) as in the simulated data
# to use the same codes, replace dataframe 'data' by 'datav2'
# datav2<- data[,c(testcols,'Recent Stat','Intermediate Stat','Distant Stat')]
# colnames(datav2)<- c(testcols, truecols); data<- datav2

# b. ROC curves
library(pROC)
rocanali <- function(esti, predi, dati) {
  auc <- matrix(0, nrow = length(predi), ncol = length(esti))
  colnames(auc) <- esti
  rownames(auc) <- predi 
  rocs <- vector("list", length(esti))
  names(rocs) <- esti  
  for (i in seq_along(esti)) {
    roc_data <- list()  
    for (j in seq_along(predi)) {
      rocobj <- roc(response = dati[[esti[i]]], predictor = dati[[predi[j]]])
      auc[j, i] <- round(rocobj$auc, 4)     
      TPR <- rocobj$sensitivities
      FPR <- 1 - rocobj$specificities     
      ff <- data.frame(FPR, TPR)     
      # Bin the data if the number of rows is too large
      if (nrow(ff) > 10000) {
        ff$bin <- cut(ff$FPR, breaks = seq(-0.0001,1,0.0001))
        ff <- aggregate(TPR ~ bin, data = ff, FUN = mean)
        ff$FPR <- (as.numeric(ff$bin) - 1) * 0.0001
        ff$bin <- NULL
      }    
      ff$pred <- predi[j]
      roc_data[[j]] <- ff 
    }  
    rocs[[i]] <- do.call(rbind, roc_data)
  }
  return(list(rocs, as.data.frame(auc)))
}

pred  <- testcols
est   <- truecols
roqui <- rocanali(est,pred,data)

# c. Extract AUC values
auc<- roqui[[2]]

# d. The proportion of true HBD SNPs in offspring as a function of locus-specific predicted HBD values
library('data.table')
classify_value <- function(value) {
  if (value >= 0 & value < 0.05) {
    return(1)
  } else if (value >= 0.2 & value < 0.3) {
    return(2)
  } else if (value >= 0.45 & value < 0.55) {
    return(3)
  } else if (value >= 0.7 & value < 0.8) {
    return(4)
  } else if (value >= 0.95 & value <= 1) {
    return(5)
  } else {
    return(NA)
  }
}
for (method in testcols) {
  data[[paste0(method, "_class")]] <- sapply(data[[method]], classify_value)
}
propTRUElist = list()
for(tester in paste0(testcols, "_class") ){
        props = as.data.frame(matrix(NA,nrow=5,ncol=5))
        props[,1] = tester
        tmp = data[,c(truecols,tester)]
        tmp = tmp[!is.na(tmp[,tester]),]
        for (k in 1:5){
           props[k,2] = k
           tmp2 = tmp[tmp[,tester] == k,]
           total_count = nrow(tmp2)
           if (total_count > 0) {
                prop_recent<- nrow(tmp2[tmp2$`Recent Base Population` == 1,])/total_count
                prop_intermediate<- nrow(tmp2[tmp2$`Intermediate Base Population` == 1,])/total_count
                prop_distant<- nrow(tmp2[tmp2$`Distant Base Population` == 1,])/total_count
                props[k,c(3,4,5)]=c(prop_recent,prop_intermediate,prop_distant)
           }
        }
        propTRUElist[[tester]] = props
        rm(props,tmp,tmp2,prop_recent,prop_intermediate,prop_distant)
}
propTrueHBD <- as.data.frame(rbindlist(propTRUElist))
propTrueHBD[,1]<- gsub("_class","",propTrueHBD[,1])
propTrueHBD[,2]<- (propTrueHBD[,2]-1)*0.25

colnames(propTrueHBD)<- c('Method','PredictedHBDclass','Recent Base Population','Intermediate Base Population','Distant Base Population')
save(corglob, corloc, meanF, roqui, auc, propTrueHBD, file='performance_sim_mediumdensity.RData')




