library(RZooRoH)

mysamples <- '../input/haplocombi_samples.txt'
frq <- read.table('../input/allele_frequencies.txt',header=FALSE)

my.data.pp <- zoodata('../input/haplocombi_pat_pat_gen.txt', samplefile = mysamples, allelefreq=frq$V1)
my.model.pp <- zoomodel(K=7, base=5, layers=TRUE)
my.resu.pp <- zoorun(my.model.pp, my.data.pp, localhbd = TRUE, nT=2)

my.data.pm <- zoodata('../input/haplocombi_pat_mat_gen.txt', samplefile = mysamples, allelefreq=frq$V1)
my.model.pm <- zoomodel(K=7, base=5, layers=TRUE)
my.resu.pm <- zoorun(my.model.pm, my.data.pm, localhbd = TRUE, nT=2)

my.data.mp <- zoodata('../input/haplocombi_mat_pat_gen.txt', samplefile = mysamples, allelefreq=frq$V1)
my.model.mp <- zoomodel(K=7, base=5, layers=TRUE)
my.resu.mp <- zoorun(my.model.mp, my.data.mp, localhbd = TRUE, nT=2)

my.data.mm <- zoodata('../input/haplocombi_mat_mat_gen.txt', samplefile = mysamples, allelefreq=frq$V1)
my.model.mm <- zoomodel(K=7, base=5, layers=TRUE)
my.resu.mm <- zoorun(my.model.mm, my.data.mm, localhbd = TRUE, nT=2)

save(my.model.pp, 
     my.data.pp, my.data.pm, my.data.mp, my.data.mm,
     my.resu.pp, my.resu.pm, my.resu.mp, my.resu.mm,
     file='ZR_Mix7L_haplocombi.RData')
quit(save='no')


