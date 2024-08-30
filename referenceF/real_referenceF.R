library(RZooRoH)
options(scipen=9999)

load('real_ZR_Mix7L.RData')

##### Compute Globals estimates
reali<- my.resu@realized 
# fields (columns) correspond to the value for each HBD class (HBD_R5, HBD_R25, HBD_R125, HBD_R625, HBD_R3126, HBD_R15625, nonHBD_R15625)
# we define three base populations: base=R25 (young), base=R125 (intermediate) and base=R3125 (distant)
# and accumulate the realized inbreeding of the HBD classes <= base

zooroh3125 <- as.data.frame(rowSums(reali[,1:5]))
zooroh125 <- as.data.frame(rowSums(reali[,1:3]))
zooroh25  <- as.data.frame(rowSums(reali[,1:2]))

true_gwF <- cbind(my.data@sample_ids,zooroh25,zooroh125,zooroh3125)
colnames(true_gwF)<- c('ind','HBDp25','HBDp125','HBDp3125')
save(true_gwF, file = 'real_glob_reference_HBDp.RData')

##### Compute Locus-specific HBD estimates (HBD probabilities)
pos <- my.data@bp
chr_vec <- vector("character", length(pos))
for (i in seq_len(nrow(my.data@chrbound))) {
    start_idx <- my.data@chrbound[i, 1]
    end_idx <- my.data@chrbound[i, 2]
    chr_vec[start_idx:end_idx] <- my.data@chrnames[i]
}
snpid = paste0(chr_vec,"_",pos)
nsnp = lenth(pos)
results_list <- list()
for ( j in seq_along(my.resu@sampleids) ){
   aniid  = my.resu@sampleids[j]
   tmp  = my.resu@hbdp[[j]]
   tmp = t(tmp)
   aniidvec       <- rep(aniid, nsnp)
   zooroh25       <- tmp[,1] + tmp[,2]
   zooroh125      <- tmp[,1] + tmp[,2] + tmp[,3]
   zooroh3125     <- tmp[,1] + tmp[,2] + tmp[,3] + tmp[,4] + tmp[,5]
   combined = cbind(aniidvec, snpid, zooroh25, zooroh125, zooroh3125)
   results_list[[length(results_list) + 1]] <- combined
}
ff <- do.call(rbind, results_list)
ff <- data.frame(ff, stringsAsFactors = FALSE)
ff[, 3] <- as.numeric(ff[, 3]); ff[, 4] <- as.numeric(ff[, 4]); ff[, 5] <- as.numeric(ff[, 5])

ftrue_HBDp<- ff
colnames(ftrue_HBDp)<- c('id','snpid','HBDp25','HBDp125','HBDp3125')
save(ftrue_HBDp, file = 'real_loc_reference_HBDp.RData')

##### Obtain estimated Locus-specific HBD class HBD class # will be used for the ROC curves
# uses results from the viterbi algorithm
results_list <- list()
vit  = my.resu@hbdseg
#dataframe with the following fields
#id, chrom, start_snp, end_snp, start_pos, end_pos, number_snp, length, HBDclass

for ( j in seq_along(my.resu@sampleids) ){
      aniid       <- my.resu@sampleids[j]
      aniidvec    <- rep(aniid, nsnp)
      stat        <- integer(nsnp)  
      class       <- integer(nsnp) 
      tmp         <- vit[vit$id == j, ]
      nseg        <- nrow(tmp)	  
      if (nseg == 0) {
           results_list[[length(results_list) + 1]] <- cbind(aniidvec, snpid, stat, class)
           next
      }
      for (k in 1:nseg){	  
              ch = my.data@chrnames[ tmp[k,'chrom'] ]
              ini  = tmp[k,'start_pos']
              fin  = tmp[k,'end_pos']
              cla  = tmp[k,'HBDclass']
	      logi <- chr_vec == ch & pos >= ini & pos <= fin
              stat[logi]  <- 1
              class[logi] <- cla
      }
      results_list[[length(results_list) + 1]] <- cbind(aniidvec, snpid, stat, class)
}

ff <- do.call(rbind, results_list)
ff <- data.frame(ff, stringsAsFactors = FALSE)
ff[, 'stat'] <- as.numeric(ff[, 'stat'])
ff[, 'class'] <- as.numeric(ff[, 'class'])

ftrue_HBDclass<- ff
colnames(ftrue_HBDclass)<- c('id','snpid','stat','class')
save(ftrue_HBDclass, file = 'real_loc_reference_HBDclass.RData')
