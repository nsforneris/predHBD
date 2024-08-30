options(scipen=9999)
library(dplyr)
library(forcats)
library(reshape2)
library(ggpubr)

# codes to extract local kinship (1 per SNP) for each couple
# files read correspond to those from the moderately inbred simulation with a medium density array
# but the codes are valid for the other data sets

# comon variables or files to more than one method ####
ped <- read.table( "../predict/input/sim_ped" )
colnames(ped) <- c('id','ids','idd')
rownames(ped) <- as.character(ped$id)
targ<- read.table( "../predict/input/sim_targets" )[,2]
targ<- ped[as.character(targ),]
rm(ped)

# snps used for prediction
nchr=25
mapfi='../predict/input/sim_nontargets.bim'
chip <- read.table(mapfi)[,c(1,2,4)]
colnames(chip)<- c('ch','snpid','pos')
rm(mapfi)

# for easy matching files
off_pairs         <- paste(targ$ids, targ$idd, sep = "_")
off_pairs_reverse <- paste(targ$idd, targ$ids, sep = "_")

#############
# predictions
#############
fpred = list()

# IBD_Haplo9c ####
coup<-read.table("../predict/ibd_haplo/output_unphased/1_unphased_2011.ids")[,c('V1','V4','V8')] #same for all chromosomes
colnames(coup)<- c('scoreset','V1','V2')
for (i in rownames(targ)){
        ps = targ[i,'ids']
        pd = targ[i,'idd']
        score= which( (coup$V1 == ps & coup$V2 == pd) | (coup$V1 == pd & coup$V2 == ps) )
        targ[i,'scoreset'] = score
}
results_list=list()
for (i in 1: nchr){
  ffdir = paste0('../predict/ibd_haplo/output_unphased/',i,'_unphased_2011.qibd')
  ff <- read.table(ffdir, header = FALSE)
  colnames(ff) <- c('scoreset','mrk','pos','S1','S2','S3','S4','S5','S6','S7','S8','S9')
  mapc <- chip[chip$ch == i, ]
  id       <- targ$id[  match(ff$scoreset,targ$scoreset) ]  
  snpid    <- mapc$snpid[ ff$mrk ]
  IBD_Haplo9c <- ff$S1 + 0.5 * (ff$S3 + ff$S5 + ff$S7) + 0.25 * ff$S8
  results_list[[i]]  <- cbind(id, snpid, IBD_Haplo9c)
}
ff <- do.call(rbind, results_list)
ff <- data.frame(ff, stringsAsFactors = FALSE)
ff[,'IBD_Haplo9c'] <- as.numeric(ff[, 'IBD_Haplo9c'])
fpred[["IBD_Haplo9c"]] = ff
#2500000 rows
rm(coup, ff, ffdir, i, IBD_Haplo9c, id, mapc, pd, ps, results_list, score, snpid)

# IBD_Haplo15c ####
coup<-read.table("../predict/ibd_haplo/output_phased/1_phased_2011.ids")[,c('V1','V4','V8')] #same for all chromosomes
colnames(coup)<- c('scoreset','V1','V2')
for (i in rownames(targ)){
        ps = targ[i,'ids']
        pd = targ[i,'idd']
        score= which( (coup$V1 == ps & coup$V2 == pd) | (coup$V1 == pd & coup$V2 == ps) )
        targ[i,'scoreset'] = score
}
results_list=list()
for (i in 1: nchr){
  ffdir = paste0('../predict/ibd_haplo/output_phased/',i,'_phased_2011.qibd')
  ff <- read.table(ffdir, header = FALSE)
  colnames(ff)<- c('scoreset','mrk','pos','S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15')
  mapc     <- chip[chip$ch == i, ]
  id       <- targ$id[  match(ff$scoreset,targ$scoreset) ]
  snpid    <- mapc$snpid[ ff$mrk ]
  s1 = ff$S1
  s3 = ff$S3  + ff$S4
  s5 = ff$S6  + ff$S7
  s7 = ff$S9  + ff$S10
  s8 = ff$S11 + ff$S12 + ff$S13 + ff$S14
  IBD_Haplo15c <- s1 + 0.5 * (s3 + s5 + s7) + 0.25*s8
  results_list[[i]]  <- cbind(id, snpid, IBD_Haplo15c)
}
ff <- do.call(rbind, results_list)
ff <- data.frame(ff, stringsAsFactors = FALSE)
ff[,'IBD_Haplo15c'] <- as.numeric(ff[, 'IBD_Haplo15c'])
fpred[["IBD_Haplo15c"]] = ff
#2500000 rows
rm(coup, ff, ffdir, i, IBD_Haplo15c, id, mapc, pd, ps, results_list, s1, s3, s5, s7, s8, score, snpid)

# GIBDLD ####
results_list <- list()
for (ch in 1:nchr) {
  spec <- read.table(paste0('../predict/GIBDLD/prefix_', ch, '.ibdtxt.gz'), nrows = 1, stringsAsFactors = FALSE)
  dim <- spec[1, 2]
  nmk <- spec[1, 4]
  spec[] <- lapply(spec, as.character) # Convert all columns to character
  loc <- read.table(paste0('../predict/GIBDLD/prefix_', ch, '.ibdtxt.gz'), skip = 1, stringsAsFactors = FALSE)
  v1_vec <- loc[, 2]
  v2_vec <- loc[, 4]
  for (j in 1:nmk) {
    pos <- 4 + j * 10
    snpid_vec <- spec[[4 + j]]
    gibdld_vec <- loc[, pos]
    combined <- cbind(v1_vec, v2_vec, snpid_vec, gibdld_vec)
    results_list[[length(results_list) + 1]] <- combined
  }
}
ff <- do.call(rbind, results_list)
ff <- data.frame(ff, stringsAsFactors = FALSE)  
ff[, 4] <- as.numeric(ff[, 4])
colnames(ff) = c('coup1','coup2','snpid','GIBDLD')
ff_pairs    <- paste(ff$coup1, ff$coup2, sep = "_")
match_idx_1 <- match(ff_pairs, off_pairs)
match_idx_2 <- match(ff_pairs, off_pairs_reverse)
ff$id[!is.na(match_idx_1)] <- targ$id[ match_idx_1[ !is.na(match_idx_1) ]]
ff$id[!is.na(match_idx_2)] <- targ$id[ match_idx_2[ !is.na(match_idx_2) ]]
ff = ff[!is.na(ff$id), c('id','snpid','GIBDLD')]
fpred[["GIBDLD"]] = ff; rm(ff)
# 2498700 rows
rm(ch, combined, dim, ff_pairs, gibdld_vec, j, loc, match_idx_1, match_idx_2, nmk, 
pos, results_list, snpid_vec, spec, v1_vec, v2_vec)

# LocalNgsRelate ####
ffdir = '../predict/LocalNgsRelate/loc_locngsrelate.txt'
ff = read.table(ffdir, header=F)
colnames(ff) = c('coup1','coup2','snpid','val','LocalNgsRelate'); ff$val=NULL
ff_pairs    <- paste(ff$coup1, ff$coup2, sep = "_")
match_idx_1 <- match(ff_pairs, off_pairs)
match_idx_2 <- match(ff_pairs, off_pairs_reverse)
ff$id[!is.na(match_idx_1)] <- targ$id[ match_idx_1[ !is.na(match_idx_1) ]]
ff$id[!is.na(match_idx_2)] <- targ$id[ match_idx_2[ !is.na(match_idx_2) ]]
ff = ff[!is.na(ff$id), c('id','snpid','LocalNgsRelate')]
fpred[["LocalNgsRelate"]] = ff
#2498700 rows
rm(ff, ffdir, ff_pairs, match_idx_1, match_idx_2)

# TRUFFLE ####
ff = '../predict/TRUFFLE/truffle.segments'
ff = read.table(ff,header=T)
colnames(ff)= c('TYPE','coup1','coup2','ch','VARSTART','VAREND','POS','Mbp','LENGTH','Mbp.1','NMARKERS')
ff_pairs    <- paste(ff$coup1, ff$coup2, sep = "_")
match_idx_1 <- match(ff_pairs, off_pairs)
match_idx_2 <- match(ff_pairs, off_pairs_reverse)
ff$id[!is.na(match_idx_1)] <- targ$id[ match_idx_1[ !is.na(match_idx_1) ]]
ff$id[!is.na(match_idx_2)] <- targ$id[ match_idx_2[ !is.na(match_idx_2) ]]
ff = ff[!is.na(ff$id), ]
ff[ff$TYPE=="IBD1",'kin']=0.25
ff[ff$TYPE=="IBD2",'kin']=0.5
nam = nrow(targ)
map = chip
for (i in 1:nam){
   map$id = targ[i,'id']
   tmp = ff[ff$id == targ[i,'id'], ]
   map$TRUFFLE=0
   tmp1 = tmp[tmp$TYPE=="IBD1",]
   tmp2 = tmp[tmp$TYPE=="IBD2",]
   nseg = nrow(tmp1)
   if (nseg > 0){
          for (k in 1:nseg){
                  ini  = tmp1[k,'VARSTART']+1; fin  = tmp1[k,'VAREND']+1
                  map[ini:fin, 'TRUFFLE']  =  tmp1[k,'kin']
          }
          rm(ini,fin)
   }
   nseg = nrow(tmp2)
   if (nseg > 0){
          for (k in 1:nseg){
                  ini  = tmp2[k,'VARSTART']+1; fin  = tmp2[k,'VAREND']+1
                  map[ini:fin, 'TRUFFLE']  =  tmp2[k,'kin']
          }
          rm(ini,fin)
   }
   if (i == 1) {
           mapprev = map
   } else {
           mapprev = rbind(mapprev,map)
   }
   rm(nseg,tmp, tmp1, tmp2)
}
mapprev$ch=NULL; mapprev$pos=NULL
fpred[["TRUFFLE"]] = mapprev
#2500000 rows
rm(ff, ff_pairs, i, k, map, mapprev, match_idx_1, match_idx_2, nam)

# ZooRoH Mix7L ####
library(RZooRoH)
load('../predict/RZooRoH/ZR_Mix7L_haplocombi.RData')
pos <- my.data.pp@bp
chr_vec <- vector("character", length(pos))
for (i in seq_len(nrow(my.data.pp@chrbound))) {
    start_idx <- my.data.pp@chrbound[i, 1]
    end_idx <- my.data.pp@chrbound[i, 2]
    chr_vec[start_idx:end_idx] <- my.data.pp@chrnames[i]
}
snpid = paste0(chr_vec,"_",pos)
match_vec = match(my.resu.pp@sampleids,off_pairs)
nam = nrow(targ)
nsnp = length(pos)
results_list <- list()
for ( j in 1:nam ){
   aniid  = targ[match_vec[j],'id']
   tmp  = ( my.resu.pp@hbdp[[j]] + my.resu.pm@hbdp[[j]] +
            my.resu.mp@hbdp[[j]] + my.resu.mm@hbdp[[j]] ) / 4
   tmp = t(tmp)
   aniidvec <- rep(aniid, nsnp)
   zooroh25      <- tmp[,1] + tmp[,2]
   zooroh125     <- tmp[,1] + tmp[,2] + tmp[,3]
   combined = cbind(aniidvec, snpid, zooroh25, zooroh125)
   results_list[[length(results_list) + 1]] <- combined
}
ff <- do.call(rbind, results_list)
ff <- data.frame(ff, stringsAsFactors = FALSE)  
ff[, 3] <- as.numeric(ff[, 3]); ff[, 4] <- as.numeric(ff[, 4])
colnames(ff)<- c('id', 'snpid', 'ZooRoH-25', 'ZooRoH-125')
fpred[["ZooRoH"]]<- ff
rm(aniid, aniidvec, chr_vec, combined, end_idx, ff, i, j, match_vec,
my.data.mm, my.data.mp, my.data.pm, my.data.pp, my.model.pp,
my.resu.mm, my.resu.mp, my.resu.pm, my.resu.pp,
nam, nsnp, pos, results_list, snpid, start_idx, tmp, zooroh125, zooroh25)

# ZooRoH 1R ####
#library(RZooRoH) 
load('../predict/RZooRoH/ZR_1R_haplocombi.RData')
pos <- my.data.pp@bp
chr_vec <- vector("character", length(pos))
for (i in seq_len(nrow(my.data.pp@chrbound))) {
    start_idx <- my.data.pp@chrbound[i, 1]
    end_idx <- my.data.pp@chrbound[i, 2]
    chr_vec[start_idx:end_idx] <- my.data.pp@chrnames[i]
}
snpid = paste0(chr_vec,"_",pos)
match_vec = match(my.resu.pp@sampleids,off_pairs)
nam = nrow(targ)
nsnp = length(pos)
results_list <- list()
for ( j in 1:nam ){
   aniid  = targ[match_vec[j],'id']
   tmp  = ( my.resu.pp@hbdp[[j]] + my.resu.pm@hbdp[[j]] +
            my.resu.mp@hbdp[[j]] + my.resu.mm@hbdp[[j]] ) / 4
   tmp = t(tmp)
   aniidvec <- rep(aniid, nsnp)
   zooroh1r      <- tmp[,1]
   combined = cbind(aniidvec, snpid, zooroh1r)
   results_list[[length(results_list) + 1]] <- combined
}
ff <- do.call(rbind, results_list)
ff <- data.frame(ff, stringsAsFactors = FALSE)  
ff[, 3] <- as.numeric(ff[, 3])
colnames(ff)<- c('id', 'snpid', 'ZooRoH-1R')
fpred[["ZooRoH-1R"]]<- ff
rm(aniid, aniidvec, chr_vec, combined, end_idx, ff, i, j, match_vec,
my.data.mm, my.data.mp, my.data.pm, my.data.pp, my.model.pp,
my.resu.mm, my.resu.mp, my.resu.pm, my.resu.pp,
nam, nsnp, pos, results_list, snpid, start_idx, tmp, zooroh1r)

# PLINK ROH ####
ff.pp=read.table(file='../predict/ROH/fake_proh_pp.hom',header=T)
ff.pm=read.table(file='../predict/ROH/fake_proh_pm.hom',header=T)
ff.mp=read.table(file='../predict/ROH/fake_proh_mp.hom',header=T)
ff.mm=read.table(file='../predict/ROH/fake_proh_mm.hom',header=T)
ff<- rbind(ff.pp, ff.pm, ff.mp, ff.mm)
ff_pairs    <- ff$IID
match_idx_1 <- match(ff_pairs, off_pairs)
ff$id[!is.na(match_idx_1)] <- targ$id[ match_idx_1[ !is.na(match_idx_1) ]]
nam = nrow(targ)
map = chip
for (i in 1:nam){
   map$id = targ[i,'id']
   tmp = ff[ff$id == targ[i,'id'], ]
   map$`PLINK ROH`=0; nseg = nrow(tmp)
   if (nseg > 0){
          for (k in 1:nseg){
                  chr  = tmp[k, 'CHR']; ini  = tmp[k,'POS1']; fin  = tmp[k,'POS2']
                  logi = map$ch == chr & map$pos>=ini & map$pos<=fin
                  map[logi, 'PLINK ROH']  = map[logi,'PLINK ROH'] + 0.25
          }
          rm(chr,ini,fin,logi)
   }
   if (i == 1) {
           mapprev = map
   } else {
           mapprev = rbind(mapprev,map)
   }
   rm(nseg,tmp)
}
mapprev$ch=NULL; mapprev$pos=NULL
fpred[["PLINK ROH"]] = mapprev
#2500000 rows
rm(ff, ff_pairs, ff.mm, ff.mp, ff.pm, ff.pp, i, k, map, mapprev, match_idx_1, nam)

# phasedibd ####
ff = '../predict/phasedibd/loc_phasedibd.txt_clean'
ff = read.table(ff)
colnames(ff)= c('pos1','pos2','coup1','coup2','hap1','hap2','ch')
ff_pairs    <- paste(ff$coup1, ff$coup2, sep = "_")
match_idx_1 <- match(ff_pairs, off_pairs)
match_idx_2 <- match(ff_pairs, off_pairs_reverse)
ff$id[!is.na(match_idx_1)] <- targ$id[ match_idx_1[ !is.na(match_idx_1) ]]
ff$id[!is.na(match_idx_2)] <- targ$id[ match_idx_2[ !is.na(match_idx_2) ]]
ff = ff[!is.na(ff$id), ]

nam = nrow(targ)
map = chip
for (i in 1:nam){
   map$id = targ[i,'id']
   tmp = ff[ff$id == targ[i,'id'], ]
   map$phasedibd=0; nseg = nrow(tmp)
   if (nseg > 0){
          for (k in 1:nseg){
                  chr  = tmp[k, 'ch']; ini  = tmp[k,'pos1']; fin  = tmp[k,'pos2']
                  logi = map$ch == chr & map$pos>=ini & map$pos<=fin
                  map[logi, 'phasedibd']  = map[logi,'phasedibd'] + 0.25
          }
          rm(chr,ini,fin,logi)
   }
   if (i == 1) {
           mapprev = map
   } else {
           mapprev = rbind(mapprev,map)
   }
   rm(nseg,tmp)
}
mapprev$ch=NULL; mapprev$pos=NULL
fpred[["phasedibd"]] = mapprev; rm(mapprev)
#2500000 rows
rm(ff, ff_pairs, i, k, map, match_idx_1, match_idx_2, nam)

# Hap-IBD ####
ff = '../predict/hap-IBD/hapibd.ibd.gz'
ff = read.table(ff)
colnames(ff)= c('coup1','hap1','coup2','hap2','ch','pos1','pos2','len')
ff_pairs    <- paste(ff$coup1, ff$coup2, sep = "_")
match_idx_1 <- match(ff_pairs, off_pairs)
match_idx_2 <- match(ff_pairs, off_pairs_reverse)
ff$id[!is.na(match_idx_1)] <- targ$id[ match_idx_1[ !is.na(match_idx_1) ]]
ff$id[!is.na(match_idx_2)] <- targ$id[ match_idx_2[ !is.na(match_idx_2) ]]
ff = ff[!is.na(ff$id), ]

nam = nrow(targ)
map = chip
for (i in 1:nam){
   map$id = targ[i,'id']
   tmp = ff[ff$id == targ[i,'id'], ]
   map$`hap-IBD`=0; nseg = nrow(tmp)
   if (nseg > 0){
          for (k in 1:nseg){
                  chr  = tmp[k, 'ch']; ini  = tmp[k,'pos1']; fin  = tmp[k,'pos2']
                  logi = map$ch == chr & map$pos>=ini & map$pos<=fin
                  map[logi, 'hap-IBD']  = map[logi,'hap-IBD'] + 0.25
          }
          rm(chr,ini,fin,logi)
   }
   if (i == 1) {
           mapprev = map
   } else {
           mapprev = rbind(mapprev,map)
   }
   rm(nseg,tmp)
}
mapprev$ch=NULL; mapprev$pos=NULL
fpred[["hap-IBD"]] = mapprev; rm(mapprev)
#2500000 rows
rm(ff, ff_pairs, i, k, map, match_idx_1, match_idx_2, nam)

# GERMLINE ####
ff = '../predict/GERMLINE/loc_germline_haplo.txt'
ff = read.table(ff)
colnames(ff)= c('fam1','coup1','hap1','fam2','coup2','hap2','ch','pos1','pos2','snpid1','snpid2','nsnp','len','unit','dat1','dat2','dat3')
ff_pairs    <- paste(ff$coup1, ff$coup2, sep = "_")
match_idx_1 <- match(ff_pairs, off_pairs)
match_idx_2 <- match(ff_pairs, off_pairs_reverse)
ff$id[!is.na(match_idx_1)] <- targ$id[ match_idx_1[ !is.na(match_idx_1) ]]
ff$id[!is.na(match_idx_2)] <- targ$id[ match_idx_2[ !is.na(match_idx_2) ]]
ff = ff[!is.na(ff$id), ]

nam = nrow(targ)
map = chip
for (i in 1:nam){
   map$id = targ[i,'id']
   tmp = ff[ff$id == targ[i,'id'], ]
   map$GERMLINE=0; nseg = nrow(tmp)
   if (nseg > 0){
          for (k in 1:nseg){
                  chr  = tmp[k, 'ch']; ini  = tmp[k,'pos1']; fin  = tmp[k,'pos2']
                  logi = map$ch == chr & map$pos>=ini & map$pos<=fin
                  map[logi, 'GERMLINE']  = map[logi,'GERMLINE'] + 0.25
          rm(chr,ini,fin,logi)
          }
   }
   if (i == 1) {
           mapprev = map
   } else {
           mapprev = rbind(mapprev,map)
   }
   rm(nseg,tmp)
}
mapprev$ch=NULL; mapprev$pos=NULL
fpred[["GERMLINE"]] = mapprev; rm(mapprev)
# 2500000 rows
rm(ff, ff_pairs, i, k, map, match_idx_1, match_idx_2, nam)

# Refined IBD ####
ff = paste('../predict/refinedibd/refinedibd.ibd.merged.gz')
ff = read.table(ff)
colnames(ff)= c('coup1','hap1','coup2','hap2','ch','pos1','pos2','lod','len')
ff_pairs    <- paste(ff$coup1, ff$coup2, sep = "_")
match_idx_1 <- match(ff_pairs, off_pairs)
match_idx_2 <- match(ff_pairs, off_pairs_reverse)
ff$id[!is.na(match_idx_1)] <- targ$id[ match_idx_1[ !is.na(match_idx_1) ]]
ff$id[!is.na(match_idx_2)] <- targ$id[ match_idx_2[ !is.na(match_idx_2) ]]
ff = ff[!is.na(ff$id), ]

nam = nrow(targ)
map = chip
for (i in 1:nam){
   map$id = targ[i,'id']
   tmp = ff[ff$id == targ[i,'id'], ]
   map$`Refined IBD`=0; nseg = nrow(tmp)
   if (nseg > 0){
          for (k in 1:nseg){
                  chr  = tmp[k, 'ch']; ini  = tmp[k,'pos1']; fin  = tmp[k,'pos2']
                  logi = map$ch == chr & map$pos>=ini & map$pos<=fin
                  map[logi, 'Refined IBD']  = map[logi,'Refined IBD'] + 0.25
          }
          rm(chr,ini,fin,logi)
   }
   if (i == 1) {
           mapprev = map
   } else {
           mapprev = rbind(mapprev,map)
   }
   rm(nseg,tmp)
}
mapprev$ch=NULL; mapprev$pos=NULL
fpred[["Refined IBD"]] = mapprev; rm(mapprev)
# 2500000 rows
rm(ff, ff_pairs, i, k, map, match_idx_1, match_idx_2, nam)

######################
# reference inbreeding
######################
ftrue=read.table("../referenceF/loc_truef.gz")
colnames(ftrue)<- c('id','snpid','node','nodeage')
ftrue[,'Recent Base Population'] = 0
ftrue[,'Intermediate Base Population'] = 0
ftrue[,'Distant Base Population'] = 0
ftrue[ftrue$nodeage<=15,'Recent Base Population'] = 1
ftrue[ftrue$nodeage<=50,'Intermediate Base Population'] = 1
ftrue[ftrue$nodeage<=500,'Distant Base Population'] = 1
ftrue[,c('node','nodeage')] = NULL

# if real data ####
  # load HBD probabilities ###
  # load("../referenceF/real_loc_reference_HBDp.RData")
  # colnames(ftrue_HBDp)<- c('id','snpid','Recent Base Population',
  #		    'Intermediate Base Population',
  #		    'Distant Population')
  # load HBD status ###
  # load("../referenceF/real_loc_reference_HBDclass.RData")
  # colnames(ftrue_HBDclass)<- c('id','snpid','stat','class')
  # ftrue_HBDclass <- ftrue_HBDclass %>%
  # mutate(
  #  `Recent Stat` = ifelse(class > 0 & class <= 2, 1, 0),
  #  `Intermediate Stat` = ifelse(class > 0 & class <= 3, 1, 0),
  #  `Distant Stat` = ifelse(class > 0 & class <= 5, 1, 0)
  # )
  # 

######################
# merge into a single table predictions and reference F
######################

merged_fpred <- Reduce(function(x, y) merge(x, y, by = c("id", "snpid")), fpred)
merged_fpred <- merged_fpred[, c('id','snpid',
                               'IBD_Haplo15c','IBD_Haplo9c','GIBDLD','LocalNgsRelate','TRUFFLE',
							   'ZooRoH-125','ZooRoH-25','ZooRoH-1R','PLINK ROH',
							   'phasedibd','hap-IBD','GERMLINE','Refined IBD')] 

data = merge(merged_fpred, ftrue, by=c('id','snpid'))

save(data, file='local_sim_mediumdensity.RData')

quit(save="no")
