## script to compute Mean and CI of the performance measure and test the significance of the difference in performance between methods
## it is ilustrated for global correlations, but the same code could be applied to other performance measures

options(scipen=9999)
library(reshape2)
library(dplyr)

truecols<- c('Recent Base Population','Intermediate Base Population',
             'Distant Base Population')
testcols<- c("IBD_Haplo15c","IBD_Haplo9c","GIBDLD","LocalNgsRelate","TRUFFLE",
             "ZooRoH-125", "ZooRoH-25","ZooRoH-1R","PLINK ROH",
             "phasedibd","hap-IBD","GERMLINE", "Refined IBD",
             "UNI","GRM","Pedigree")

# integrate results from all replicates 
load('corglob_sim_mediumdensity_allreplicates.RData')

fold <- paste('r', seq(1, 100), sep="")
results_list <- vector("list", length = 100)
for (rr in seq_along(fold)) {
  ddc <- as.data.frame( corglob[[rr]] )[testcols, truecols]
  ddc$modc <- testcols
  ddc.melt <- melt(ddc, id.vars = "Method", measure.vars = truecols,
                   variable.name = "ReferenceF", value.name = "y") # as "y" could be an AUC or correlation value
  ddc.melt$rep <- fold[rr]
  results_list[[rr]] <- ddc.melt
}
dfall <- do.call(rbind, results_list)
dfall[] <- lapply(dfall, function(x) if(is.factor(x)) as.character(x) else x)
rm(fold, rr, ddc, ddc.melt, results_list)

# Mean and confidence intervals (CI) for the performance of each method  
# functions and code bellow adapted from M. Schrauf's https://github.com/schrauf/AccuracyComparer/
source("functions.R")

# fields: ReferenceF/ Method/ Mean performance/ lower CI/ upper CI/ letter given after testing the significance of the difference in performance between methods
data1 <- dfall |>
  group_by(ReferenceF) |>
  group_modify(ff)

# Test the significance of the difference in performance between pairs of methods
# fields: "ReferenceF" "Pairwise Comparison" "Mean of the difference in performance" "lower CI" "upper CI" "significance"(TRUE/FALSE) 
data2 <- dfall |>
  group_by(ReferenceF) |>
  (\(df) inner_join(df,df, by=c("ReferenceF","rep")))() |>
  filter(modc.x < modc.y) |>
  mutate(diff = y.y-y.x, comp=paste(modc.y,modc.x,sep=" - ")) |>
  group_by(ReferenceF, comp) |>
  summarise(boot_fun(diff))|>
  mutate(signif = (0 < lower) | (0 > upper))

save(data1, data2, file='corglob_sim_mediumdensity_pairdiff.RData')

