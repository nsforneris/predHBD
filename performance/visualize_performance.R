# script to plot performance results
# for one replicate of the simulation (same plots are valid for real data)

options(scipen=9999)
library(reshape2)
library(ggplot2)
library(cowplot)

load('performance_sim_mediumdensity.RData')

colo<- c("IBD_Haplo15c"="#b775e1","IBD_Haplo9c"="#8c0da4","GIBDLD"="#f8cca6","LocalNgsRelate"="#f17a74","TRUFFLE"="#efd453",
  "ZooRoH-125"="#2f4285","ZooRoH-25"="#2282f5","ZooRoH-1R"="#20d8fd","PLINK ROH"="#4f8522",
  "phasedibd"="#ec102f","hap-IBD"="#b70d61","GERMLINE"="#fd8f20","Refined IBD"="#683d0d",
  "GRM"="#938073","UNI"="#d5c5a1","Pedigree"="#d4d4d4")
truecols<- c('Recent Base Population','Intermediate Base Population',
	     'Distant Base Population')
testcols<- c("IBD_Haplo15c","IBD_Haplo9c","GIBDLD","LocalNgsRelate","TRUFFLE",
  "ZooRoH-125", "ZooRoH-25","ZooRoH-1R","PLINK ROH",
  "phasedibd","hap-IBD","GERMLINE", "Refined IBD",
  "UNI","GRM","Pedigree") 

# Correlations between predicted and reference genome-wide inbreeding
df = corglob
df = as.data.frame(df)
df$Method = rownames(df)
df.melt   <- melt(df, id.vars = c("Method"), measure.vars = truecols,
                   variable.name = "ReferenceF", value.name = "Corr")
df.melt$Method <- factor(df.melt$Method, levels = testcols)
df.melt$ReferenceF <- factor(df.melt$ReferenceF, levels = rev(truecols))

p_globcor<- ggplot(df.melt, aes(y=Corr,x=Method,color=Method,fill=Method)) + 
	    #ggplot(df.melt, aes(ymin=lower, y=Corr, ymax=upper, x=Method, color=Method, fill=Method )) +
  geom_bar(stat="identity", position = position_dodge(width = 0.9)) +
  #geom_errorbar(aes(ymin=lower,ymax=upper), position = position_dodge(width=0.9),width=0.4,size=0.6,alpha=0.9, colour="grey30") +
  facet_grid(~ReferenceF) +
  labs(title="Genome-wide Correlations", 
       subtitle="Moderate Inbreeding Scenario - Medium Density", y = "Correlation", x = " ") +
  scale_fill_manual(values = colo) + scale_color_manual(values = colo) +
  coord_cartesian(ylim = c(0.4,0.85))+
  #scale_y_continuous(breaks=seq(0.4,0.85,0.05)) + 
  theme(
    plot.title = element_text(size = 13, hjust = 0.5) , 
    plot.subtitle = element_text(size = 12, hjust = 0.5) , 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.background = element_rect(fill=NA),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(),
    panel.spacing.x = unit(15, "pt", data=NULL),
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 12, vjust = 1.5) , # grid subtitle
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.line.y = element_line(colour = "black"),
    legend.position = "none")

color_plot<- ggplot(df.melt, aes(y = Corr, x = Method, fill = Method )) +
  geom_bar(color="white",stat="identity",position=position_dodge(width = 0.9)) +
  scale_fill_manual(values = colo) +
  guides(fill = guide_legend(title="",reverse=FALSE, ncol=8)) +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.text = element_text(size=11, margin = margin(r = 5, unit = "pt")))
color_legend <- cowplot::get_legend(color_plot)

# Average predicted and reference genome-wide inbreeding levels
ef<- as.data.frame(meanF)
ef$Method<- rownames(ef)
ef$Method<- factor(ef$Method, levels = testcols)

p_meanF<- ggplot(ef[testcols,], aes(y = meanF, x = Method, fill = Method)) +
  geom_bar(stat="identity", color = "white", position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept=ef['Recent Base Population','meanF']),colour="black",lty="solid") +
  geom_hline(aes(yintercept=ef['Intermediate Base Population','meanF']),colour="black",lty="dashed") +
  geom_hline(aes(yintercept=ef['Distant Base Population','meanF']),colour="black",lty = "dotted") +
  annotate("text", x=Inf, y=ef['Recent Base Population','meanF'],label="Recent Base Population", hjust = 1.1, vjust = -0.5, size = 3.5) +
  annotate("text", x=Inf, y=ef['Intermediate Base Population','meanF'],label="Intermediate Base Population", hjust = 1.1, vjust = -0.5, size = 3.5) +
  annotate("text", x=Inf, y=ef['Distant Base Population','meanF'],label="Distant Base Population",hjust = 1.1, vjust = -0.5, size = 3.5) +
  labs(title="Average Inbreeding Levels", subtitle="Moderate Inbreeding Scenario - Medium Density",
       y="Mean Inbreeding Prediction", x = " ") +
  scale_fill_manual(values = colo) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 0.25), linetype = 'solid') +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 0.25), breaks = seq(0, 0.25, 0.05)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.ticks.y = element_line(),
    legend.position = "none"  )

# Correlations between predicted and reference locus-specific inbreeding
ff = corloc
ff = as.data.frame(ff)
ff$Method = rownames(ff)
ff.melt   <- melt(ff, id.vars = c("Method"), measure.vars = truecols,
                   variable.name = "ReferenceF", value.name = "Corr")
ff.melt$Method <- factor(ff.melt$Method, levels = testcols)
ff.melt$ReferenceF <- factor(ff.melt$ReferenceF, levels = rev(truecols))

p_loccor<- ggplot(ff.melt, aes(y=Corr,x=Method,color=Method,fill=Method)) +
            #ggplot(ff.melt, aes(ymin=lower, y=Corr, ymax=upper, x=Method, color=Method, fill=Method )) +
  geom_bar(stat="identity", position = position_dodge(width = 0.9)) +
  #geom_errorbar(aes(ymin=lower,ymax=upper), position = position_dodge(width=0.9),width=0.4,size=0.6,alpha=0.9, colour="grey30") +
  facet_grid(~ReferenceF) +
  labs(title="Locus-specific Correlations",
       subtitle="Moderate Inbreeding Scenario - Medium Density", y = "Correlation", x = " ") +
  scale_fill_manual(values = colo) + scale_color_manual(values = colo) +
  coord_cartesian(ylim = c(0.01,0.5))+
  #scale_y_continuous(breaks=seq(0,0.5,0.05)) +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5) ,
    plot.subtitle = element_text(size = 12, hjust = 0.5) ,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.background = element_rect(fill=NA),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(),
    panel.spacing.x = unit(15, "pt", data=NULL),
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 12, vjust = 1.5) , # grid subtitle
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.line.y = element_line(colour = "black"),
    legend.position = "none")

# AUC values
gf = auc
gf = as.data.frame(gf)
gf$Method = rownames(gf)
gf.melt   <- melt(gf, id.vars = c("Method"), measure.vars = truecols,
                   variable.name = "ReferenceF", value.name = "AUC")
gf.melt$Method <- factor(gf.melt$Method, levels = testcols)
gf.melt$ReferenceF <- factor(gf.melt$ReferenceF, levels = rev(truecols))

p_auc<- ggplot(gf.melt, aes(y = AUC, x = Method, color = Method )) +
  geom_point(stat="identity", size = 1.5 ) +
  facet_grid(~ ReferenceF) +
  #geom_errorbar(aes(ymin=lower, ymax=upper, color=Method), position=position_dodge(width = 0.9),
  #              width=0.8, size=0.6, alpha=0.9) +
  coord_cartesian(ylim = c(0.5,0.9))+
  scale_y_continuous(breaks=seq(0.5,0.9,0.05)) + 
  scale_color_manual(values = colo) +
  labs(title="Area Under The ROC Curves",
       subtitle="Moderate Inbreeding Scenario - Medium Density", y = "AUC", x = " ") +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5) ,
    plot.subtitle = element_text(size = 12, hjust = 0.5) ,
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.background = element_rect(fill=NA),
    panel.border = element_rect(colour = "black", fill = NA),
    #panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.spacing.x = unit(15, "pt", data=NULL),
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 12, vjust = 1.5) , # grid subtitle
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.line.y = element_line(colour = "black"),
    legend.position = "none")

# ROC curves
theme_set(theme_minimal())

rocplots = list()
for (i in truecols){
  ff = roqui[[1]][[i]]
  ff = ff[ff$pred %in% testcols,]
  auctable = gf.melt[gf.melt$ReferenceF==i,] #see previous AUC plot section
  rownames(auctable)<- auctable$Method
  # order methods according to their AUC
  ord = order(auctable$AUC, decreasing=T)
  temp <- auctable[ rownames(auctable)[ord],]
  # plot roc curves and highlight the methods with best correlation values
  colo_ord<- colo
  hl <- c("ZooRoH-125", "ZooRoH-25", "ZooRoH-1R", "IBD_Haplo9c", "IBD_Haplo15c", "PLINK ROH", "phasedibd", "hap-IBD")
  zr <- rownames(temp)[rownames(temp) %in% hl]
  nzr <- rownames(temp)[!(rownames(temp) %in% hl)]
  colo_ord[nzr]<- "grey80"
  q <- ggplot(ff) +
    geom_line( aes(FPR, TPR, group = pred, color = pred), size = 0.9 ) +  
    scale_colour_manual(values = colo_ord, breaks = zr, labels = paste0(seq(1,8), ": ", zr)) +
    labs(title="ROC\nModerate Inbreeding Scenario - Medium Density", 
	 subtitle = i,
	   x = "False Positive Rate", y = "True Positive Rate") +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5) ,
      plot.subtitle = element_text(size = 11.5, hjust = 0.5) ,
      plot.margin = margin(t = 5, b = 5, r = 5, l = 5),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      panel.grid = element_blank(),
      axis.text = element_text(size = 11),
      axis.ticks = element_line(),
      legend.text = element_text(size = 9),
      legend.title = element_blank(),
      legend.key.height = unit(0.8, 'lines'),
      legend.position = c(0.72, 0.35), # adjust as needed
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.key = element_rect(fill = NA, color = NA),
      legend.background = element_rect(color = NA, fill = "white"))
   rocplots[[i]]<- q
}
rm(ff, auctable, ord, temp, q, zr, nzr, hl, i, colo_ord)

p_rocs1<- rocplots[['Distant Base Population']] + 
         theme(plot.title=element_text(color="transparent", size=12, hjust = 0.5),
               axis.title.y = element_text(size = 11, vjust = 1.5)) 
p_rocs2<- rocplots[['Intermediate Base Population']] +
         theme(plot.title=element_text(size=12, hjust = 0.5),
               axis.title.y = element_text(color="transparent", size = 11, vjust = 1.5))  
p_rocs3<- rocplots[['Recent Base Population']] + 
         theme(plot.title=element_text(color="transparent", size=12, hjust = 0.5),
               axis.title.y = element_text(color="transparent",size = 11, vjust = 1.5)) 
p_rocs<- cowplot :: plot_grid(p_rocs1, p_rocs2, p_rocs3, nrow=1, ncol = 3)  

# Ranking

# ranking for global correlations
mtds<- testcols; rglob<- t(corglob); tt<- rglob
for (i in 1:nrow(tt)){
  values = as.numeric(rglob[i, mtds ])
  tt[i, mtds] = rank( -values, na.last = T )
}
rglob_r<- tt
rglob_r_melted <- melt(rglob_r, variable.name="Method", value.name="Rank")
colnames(rglob_r_melted)<- c('ReferenceF','Method','Rank')
rank_glob<- merge(df.melt, rglob_r_melted, by=c('Method','ReferenceF'))
rank_glob$meausure<- "Genome-wide Correlations"
rm(rglob,tt,mtds,rglob_r_melted)

# ranking for locus-specific correlations
mtds<- testcols[!(testcols %in% c('GRM','UNI','Pedigree'))]; rloc<- t(corloc[mtds,]); tt<- rloc
for (i in 1:nrow(tt)){
  values = as.numeric(rloc[i, mtds ])
  tt[i, mtds] = rank( -values, na.last = T )
}
rloc_r<- tt
rloc_r_melted <- melt(rloc_r, variable.name="Method", value.name="Rank")
colnames(rloc_r_melted)<- c('ReferenceF','Method','Rank')
rank_loc<- merge(ff.melt, rloc_r_melted, by=c('Method','ReferenceF'))
rank_loc$meausure<- "Locus-specific Correlations"
rm(rloc,tt,rloc_r_melted)

# we can do it also for AUC values
# mtds<- testcols[!(testcols %in% c('GRM','UNI','Pedigree'))]; rauc<- t(auc[mtds,]); tt<- rauc
# for (i in 1:nrow(tt)){
#  values = as.numeric(rauc[i, mtds ])
#  tt[i, mtds] = rank( -values, na.last = T )
# }
# rauc_r<- tt
# rauc_r_melted <- melt(rauc_r, variable.name="Method", value.name="Rank")
# colnames(rauc_r_melted)<- c('ReferenceF','Method','Rank')
# rank_auc<- merge(gf.melt, rauc_r_melted, by=c('Method','ReferenceF'))
# rank_auc$meausure<- "AUC values"
# rm(rauc,tt,mtds,rauc_r_melted)

join_df <- rbind(rank_glob,rank_loc)
# to show correlations values with two decimal digits and x10-1
format_sci <- function(x) {
  formatted <- sprintf("%.2f", 10*x)
  return(formatted)
}
join_df$Corr1<- format_sci(join_df$Corr)

p_rank1 <- ggplot(join_df[join_df$meausure=="Genome-wide Correlations", ],
		  aes(x = Rank, y = ReferenceF, fill = Method, color = Method)) +
  geom_tile(width = 0.9, height = 0.9, size = 0.25) +  
  geom_text(aes(label = ifelse(Method %in% c("GIBDLD","UNI","Pedigree","TRUFFLE","ZooRoH-1R"), 
                               Corr1,"")), vjust = 0.5, hjust = 0.5, size = 4, color = "grey20") +
  geom_text(aes(label = ifelse(!(Method %in% c("GIBDLD","UNI","Pedigree","TRUFFLE","ZooRoH-1R")), 
                               Corr1,"")), vjust = 0.5, hjust = 0.5, size = 4, color = "white") +
  scale_fill_manual(values = colo) + scale_color_manual(values = colo) +
  labs(title = "Genome-wide Correlation Ranking", 
       subtitle = "Moderate Inbreeding Scenario - Medium Density",
  y = "Reference Inbreeding", x = "Ranking", fill = "Method") +
  scale_x_continuous(breaks = 1:16) + coord_cartesian(expand = FALSE) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 12, hjust = 0.5) ,
        plot.subtitle = element_text(size = 11, hjust = 0.5) ,
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
 		axis.title.x = element_text(size = 10),
		axis.title.y = element_text(color = "transparent", size = 10),
		axis.text = element_text(size = 10),
        axis.ticks.y = element_blank(), axis.ticks.x = element_line(),
        legend.position = "none")

p_rank2 <- ggplot(join_df[join_df$meausure=="Locus-specific Correlations", ],
		  aes(x = Rank, y = ReferenceF, fill = Method, color = Method)) +
  geom_tile(width = 0.9, height = 0.9, size = 0.25) +  
  geom_text(aes(label = ifelse(Method %in% c("GIBDLD","UNI","Pedigree","TRUFFLE","ZooRoH-1R"), 
                               Corr1,"")), vjust = 0.5, hjust = 0.5, size = 4, color = "grey20") +
  geom_text(aes(label = ifelse(!(Method %in% c("GIBDLD","UNI","Pedigree","TRUFFLE","ZooRoH-1R")), 
                               Corr1,"")), vjust = 0.5, hjust = 0.5, size = 4, color = "white") +
  scale_fill_manual(values = colo) + scale_color_manual(values = colo) +
  labs(title = "Locus-specific Correlation Ranking",
       subtitle = "Moderate Inbreeding Scenario - Medium Density",
  y = "Reference Inbreeding", x = "Ranking", fill = "Method") +
  scale_x_continuous(breaks = 1:16) + coord_cartesian(expand = FALSE) + 
  theme_minimal() +
  theme(plot.title = element_text(size = 12, hjust = 0.5) ,
        plot.subtitle = element_text(size = 11, hjust = 0.5) ,
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
		axis.title.x = element_text(size = 10),
		axis.title.y = element_text(color = "transparent", size = 10),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.ticks.y = element_blank(), axis.ticks.x = element_line(),
        legend.position = "none")
		
p_rank<- plot_grid(p_rank1, p_rank2, rel_widths =  c(3.6,2.2)) 
p_rank_w_leg <- plot_grid(p_rank, color_legend, ncol = 1, rel_heights = c(0.95, 0.05))


# Proportion of true HBD SNPs in offspring as a function of locus-specific predicted HBD values


