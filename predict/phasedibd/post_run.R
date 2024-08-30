# cleans final output by searching for reported overlapping segments

options(scipen=9999)
library(data.table)
library(stringr)

segs <- fread('loc_phasedibd.txt')
segs<- segs[,c(1,2,3,4,11,9,10)]
segs$combi <- paste(segs$V1, segs$V2, segs$V3, segs$V4, segs$V11, sep="@")
names(segs)[names(segs)=="V9"] <- "start"
names(segs)[names(segs)=="V10"] <- "end"
segs<- segs[,c('start','end','combi')]

# repeated combinations
dupli<- unique(segs[duplicated(segs$combi),'combi'])
dupli<- dupli$combi

# one-time-represented combinations
segs_uni<- segs[!(segs$combi %in% dupli),]
names(segs_uni)[names(segs_uni)=="combi"] <- "combin"

# merge overlapping segments in repeated combinations
segs_dupli <- segs[segs$combi %in% dupli,]
rm(segs)
gc()

dim = length(dupli)
dfs<- vector("list", length = dim)
for (i in 1:dim){
  combin = dupli[i]
  testDT = segs_dupli[combi==combin,]
  df<- testDT[order(start)
  ][, .(start=min(start), end=max(end)),
    by=.(group=cumsum(c(1, tail(start, -1) > head(end, -1))))]
  dfs[[i]]<- cbind(combin, df)
}
df = rbindlist(dfs)
print('overlapping segments merged')
rm(segs_dupli, dfs, combin, dupli, testDT)
gc()

# joining unique and merged segments into a single table
df<- rbind(df[,c('start','end','combin')],segs_uni)
split_data <- str_split_fixed(df$combin, "@", 5)
df <- cbind(df, split_data)
df$combin = NULL
rm(segs_uni, split_data)

# print clean segments
write.table(df,file='loc_phasedibd.txt_clean',col.names=F,row.names=F,quote=F)

