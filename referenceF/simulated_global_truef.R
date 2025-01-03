# R program for computing the true proportiong of HBD 
# defined based on three base populations 

options(scipen=9999)
library(reshape2)
library(dplyr)
library(data.table)

glen = 2500 # in MB

age <- fread('../predict/input/out_nodeage.txt.gz')
setnames(age, c('node', 'age'))
age[, age := round(age)]

segments <- fread('../predict/input/out_segments.txt.gz', 
                   colClasses = c('integer', 'integer', 'integer', 'numeric', 'numeric', 'integer'))
setnames(segments, c('ind', 'nod1', 'nod2', 'ini', 'fin', 'node'))
segments[, len := (fin - ini) / 1e6]
segments[, `:=`(ini = NULL, fin = NULL)]

segs <- segments[, .(len = sum(len)), by = .(ind, node)]
segs_age<- left_join(segs, age, by ='node')
segs <- segs_age[, .(len = sum(len)), by = .(ind, age)]
segs[, referenceF := fifelse(age <= 15, "15G",
                             fifelse(age <= 50, "50G",
                                     fifelse(age <= 500, "500G", "nonHBD")))]

segs1 <- segs[, .(len = sum(len)), by = .(ind, referenceF)]
segs1[, len := len / glen]
segs1[,ind := paste0("tsk_",ind)]
rm(segs, segments)

# accumulate HBD over 15, 50 and 500 generations
library(reshape)
tt = cast(segs1, ind~referenceF, value='len')
tt[is.na(tt)]=0
tt[,'50G'] = tt[,'15G'] + tt[,'50G']
tt[,'500G'] = tt[,'50G'] + tt[,'500G']
tt[,c('nonHBD')] = NULL
rm(segs1)

true_gwF = tt[,c('ind','15G','50G','500G')]
save(true_gwF, file=paste0("glob_truef.RData"))

quit(save="no")

