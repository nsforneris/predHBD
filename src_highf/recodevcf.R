options(scipen=9999)
library(data.table)

# Read the input VCF file
vcf <- fread("output.vcf", header = TRUE)

# Find the line number where #CHROM is located
rr <- which(vcf$#CHROM == "#CHROM")

# Split the data into two parts
f0 <- vcf[1:(rr - 1), ]
f1 <- vcf[rr:length(vcf), ]

# Recode chromosome and position fields
breaks <- seq(0, 2500000000, by = 100000000)
f1$chr <- cut(f1$POS, breaks = breaks, labels = 1:25, include.lowest = TRUE)
f1$posr <- (f1$POS - 1) %/% 100000000 + 1

f1$#CHROM <- f1$chr
f1$POS <- f1$posr

f1[,c('chr','posr')] = NULL

# Write the modified data to a new VCF file
write.table(f1, file = "genotypes.vcf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

quit(save="no")


