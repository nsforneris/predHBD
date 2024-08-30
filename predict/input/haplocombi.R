library('data.table')
library('readr')

args= commandArgs(trailingOnly=T)

vcf_file  <- args[1]
ped_file  <- args[2]

# Read input
data_vcf <- fread(vcf_file, header = T)
genos <- as.data.frame(data_vcf[, 10:ncol(data_vcf)])

# read offspring and its parents
pairs <- read.table(ped_file); pairs$V1=NULL
pairs$pairname<- paste0(pairs[,1],"_",pairs[,2])
n_pairs<- nrow(pairs)
write.table(file = 'haplocombi_samples.txt', pairs$pairname,
         row.names=F, col.names=F, quote=F)

# map_info
map_info <- as.data.frame(data_vcf[, c(1, 3, 2, 4, 5)])
map_info[,1] = gsub("chr","",map_info[,1])
map_info[,2]<- paste0(map_info[,1],"_",map_info[,3])
rm(data_vcf)

# Split phased genotypes
left_side <- lapply(genos, function(col) sapply(strsplit(as.character(col), "\\|"), "[[", 1))
right_side <- lapply(genos, function(col) sapply(strsplit(as.character(col), "\\|"), "[[", 2))
rm(genos)

# Create data frames from split
paternal_hap <- as.data.frame(lapply(left_side, as.integer))
maternal_hap <- as.data.frame(lapply(right_side, as.integer))
names(paternal_hap) <- gsub("^X", "", names(paternal_hap))
names(maternal_hap) <- gsub("^X", "", names(maternal_hap))
rm(left_side, right_side)

# Recode to count reference alleles as 0/1/2
paternal_hap <- 1 - paternal_hap
maternal_hap <- 1 - maternal_hap
gc()

for (combo in c("pat_pat", "mat_mat", "pat_mat", "mat_pat")) {
  result <- data.frame(matrix(0, nrow = nrow(paternal_hap), ncol = n_pairs))
  if(combo == "pat_pat"){
    sire_gamete <- paternal_hap   #fake calf inherits paternal allele from sire
    dam_gamete <- paternal_hap    #fake calf inherits paternal allele from dam
  } else if (combo == "mat_mat") {
    sire_gamete <- maternal_hap
    dam_gamete <- maternal_hap
  } else if (combo == "pat_mat") {
    sire_gamete <- paternal_hap
    dam_gamete <- maternal_hap
  } else if (combo == "mat_pat") {
    sire_gamete <- maternal_hap
    dam_gamete <- paternal_hap
  }
  for (i in 1:n_pairs) {
    sire_col_name <- pairs[i,1]
    dam_col_name <- pairs[i,2]
    pat_col <- sire_gamete[[sire_col_name]]
    mat_col <- dam_gamete[[dam_col_name]]
    result[, i] <- pat_col + mat_col
  }
  result = cbind(map_info, result)
  gc()

  fwrite(file = paste0("haplocombi_",combo,"_gen.txt"), result,
         row.names=F, col.names=F, quote=F, sep = " ")
  rm(sire_gamete, dam_gamete, result)

}

quit(save="no")
