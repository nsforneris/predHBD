library(RZooRoH)

# Rscript ZR_allelefreq.R '../input/sim_nontargets_gen.txt' '../input/sim_nontargets_samples.txt' '../input/allele_frequencies.txt'
args= commandArgs(trailingOnly=T)

mygtfile  <- args[1]
mysamples <- args[2]
outprefix <- args[3]

my.data <- zoodata(mygtfile, samplefile = mysamples)
write.table(my.data@freqs, file=outprefix, col.names=F, row.names=F, quote=F)
quit(save='no')

