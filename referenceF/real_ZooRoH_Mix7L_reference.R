# script to run RZooRoH using real data to estimate reference inbreeding
# the same script was used for prediction (ZR_Mix7L_haplocombi.R)
# except that here we run RZooRoH only once for the real set of individuals 
# instead of running it 4 times (one time for each combination of parental haplotypes

library(RZooRoH)
mygtfile <- "seq_gen.txt"      # gen file for the real set of animals
mysamples <- "samples.txt"     # list with the names(ID) of the samples

my.data <- zoodata(mygtfile, samplefile = mysamples)
my.model <- zoomodel(K=7, base=5, layers=TRUE)

my.resu <- zoorun(my.model, my.data, localhbd = TRUE, nT=2)  
save(my.resu,file=resu)

# or you could do several runs to partition results into smaller objects (for instance, when using with sequence data)
# my.resu1 <- zoorun(my.model, my.data, ids= list1,localhbd = TRUE, nT=2)  
# my.resu2 <- zoorun(my.model, my.data, ids= list2,localhbd = TRUE, nT=2)  
# etc
# where list1 and list2 are the positions in the samples.txt file of desired individuals

save(my.data, my.resu,file='real_ZR_Mix7L.RData')

#save(my.data, my.resu1,file='real_ZR_Mix7L_1.RData')
#save(my.data, my.resu2,file='real_ZR_Mix7L_2.RData')

quit(save='no')
