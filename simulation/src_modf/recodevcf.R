options(scipen=9999)
library(data.table)

# Read the input VCF file
vcf <- fread("output.vcf", header = TRUE)

# Recode chromosome and position fields
positions <- seq(0, 2499999999, by = 100000000)
chromosomes <- seq(1,length(positions))
chr<- vcf$`#CHROM`
posr<- vcf$`POS`
posr2<- posr
for (i in 1:length(positions)) {
  chr[posr >= positions[i] & posr < positions[i] + 100000000 - 1] <- chromosomes[i]
  posr2[posr >= positions[i] & posr < positions[i] + 100000000 - 1] <- posr[posr >= positions[i] & posr < positions[i] + 100000000 - 1] - positions[i] + 1
}
vcf$`#CHROM` <- chr
vcf$`POS` <- posr2; rm(posr, posr2)
vcf$ID = paste0("var",vcf$ID)

# remove polymorphic sites (> 2 alleles)
vcf<- vcf[!grepl(",",vcf$ALT),]

# Write the modified data to a new VCF file
write.table(vcf[,c('ID','REF','ALT')], file = "genotypes_ref_alt", quote = FALSE, row.names = FALSE, col.names = TRUE)
line <- as.integer(system("awk '$1==\"#CHROM\"{print NR - 1}' output.vcf", intern = TRUE))
write.table(vcf, file = "genotypes.vcf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

command <- paste0("head -n ", line, " output.vcf | awk '$1!~/contig/' > header")
system(command)
command <- paste0("cat header genotypes.vcf > genotypes_tidy.vcf")
system(command)

quit(save="no")
