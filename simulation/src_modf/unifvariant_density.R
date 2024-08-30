options(scipen=9999)

args = commandArgs(trailingOnly=T)

set.seed(args[1])

df<- read.table('freq.frq',header=T)
library(dplyr)
dfbin<- df %>% mutate(MAFbin = cut(MAF, breaks=c(0, 0.05, 0.10, 0.15,
0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1)))

# uniform maf distribution
nmarkers<-min(table(dfbin[,'MAFbin']))
binclass = unique(dfbin$MAFbin)
unisnp<- list()
for (i in binclass){
        unisnp[[i]] <- sample(dfbin[dfbin$MAFbin == i,'SNP'],
        nmarkers, replace = FALSE, prob = NULL)
        write.table(file='out_unifvariants.txt', as.data.frame(unisnp[[i]]), row.names=F, col.names=F, quote=F, append=T)
}

rm(unisnp, df)

# thinning 
varunif<- read.table('out_unifvariants.txt')
colnames(varunif)<- c('id')

varinfo<- read.table('plink2.bim')
colnames(varinfo)<- c('ch','id','val','pos','a1','a2')

varunif<- merge(varunif,varinfo, by="id")
varinfo<- varinfo[,c('id','ch','val','pos','a1','a2')]

for (j in 1:25){
	bigvars0 = varinfo[varinfo$ch == j ,]
	vars0 = varunif[varunif$ch == j ,]
	for(i in 1:1000){
		ini= 100000*(as.numeric(i)-1)+1
		end= 100000*as.numeric(i)
		bigvars = bigvars0[bigvars0$pos>= ini & bigvars0$pos<= end,]
		vars = vars0[vars0$pos>= ini & vars0$pos<= end,]
                
		midi = (ini + end)/2
	        
		if(nrow(vars)>0){
			mini = min(abs(vars$pos-midi))
			snp = which(abs(vars$pos-midi)==mini)
	                write.table(file='out_unifvariants_10s_mb.txt', vars[snp[1],'id'], row.names=F, col.names=F, quote=F, append=T)
		} else if (nrow(bigvars)>0) {
			mini = min(abs(bigvars$pos-midi))
                        snp = which(abs(bigvars$pos-midi)==mini)
	                write.table(file='out_unifvariants_10s_mb.txt', bigvars[snp[1],'id'], row.names=F, col.names=F, quote=F, append=T)
		} else {
			write.table(file='missing.txt', cbind(j,i), row.names=F, col.names=F, quote=F, append=T)
		}
	}
}


rm(mini, snp, midi, ini, end, bigvars, vars, bigvars0, vars0, varinfo)

missing<- read.table('missing.txt')
vars<- read.table('out_unifvariants_10s_mb.txt')

chrs  = unique(missing$V1)
notused = varunif[!(varunif$id %in% vars$V1),]

for (i in chrs){
	nmiss = nrow(missing[missing$V1==i,])
        ids = sample(notused[notused$ch==i, 'id'], nmiss, replace = FALSE, prob = NULL)
        write.table(file='out_unifvariants_10s_mb.txt', ids, row.names=F, col.names=F, quote=F, append=T)
}

quit(save='no')

