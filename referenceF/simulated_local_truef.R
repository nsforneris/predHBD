# R program for computing the true HBD status at each loci
# defined based on three base populations

options(scipen=9999)
library(reshape2)
library(dplyr)
library(data.table)

age <- fread('../predict/input/out_nodeage.txt.gz')
setnames(age, c('node', 'age'))
age[, age := round(age)]
seg <- fread('../predict/input/out_segments.txt.gz',
                   colClasses = c('integer', 'integer', 'integer', 'numeric', 'numeric', 'integer'))
setnames(seg, c('ind', 'nod1', 'nod2', 'ini', 'fin', 'node'))
seg[, `:=`(nod1 = NULL, nod2 = NULL)]
# positions in seg are coded for a single big chromosome
# they were generated before running recode_vcf.R

outfi<- "loc_truef"

# variant info #
system("zcat ../predict/input/out_variant_info.txt.gz | awk '{print $1, $2}' > vv")
varinfo<- fread('vv', colClasses = c('numeric', 'integer'))
setnames(varinfo, c('pos','ord')) # pos matches that of seg table
varinfo[, name := paste0('var', ord)]
varinfo[, `:=`(ord = NULL)]

# specific variants that will be used for predictions #
vcf<- fread('../predict/input/plink.vcf.gz') [,1:3]
vcf[, snpid := paste0(`#CHROM`, "_", POS)]
vcf[, `:=`(`#CHROM` = NULL, POS = NULL)] 
varinfo = merge(varinfo,vcf,by.x = 'name', by.y = 'ID')
rm(vcf)

nsnp  <- nrow(varinfo)
posi  <- varinfo$pos
snpid <- varinfo$snpid
results <- data.table()

for (j in 0:99){
        seg0 <- seg[ind==j]
        ini  <- seg0$ini
        fini <- seg0$fin
        node <- seg0$node
        resi <- sapply(posi, function(pos) node[which(ini <= pos & fini > pos)])
        resiage <- age[resi + 1,]
        ani <- rep(j, nsnp)
        result <- data.table(ani, snpid, resiage)
        results <- rbind(results, result)
}

results[, ani := paste0('tsk_',ani)]
fwrite(results, outfi, row.names=F, col.names=F, quote = FALSE)
quit(save="no")



