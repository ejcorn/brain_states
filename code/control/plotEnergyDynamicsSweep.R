# relate group level control energy metrics to empirically obtained brain state dynamics
# namely, we expect persistence energy to explain persistence probability and dwell time

args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

masterdir <- paste(basedir,'results/',name_root,'/',sep='')
savedir <- paste(masterdir,'analyses/control_energy/',sep='')

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))
source(paste(basedir,'code/miscfxns/statfxns.R',sep=''))

scanlab <- c("RestComb","nBackComb")
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  
clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")

# compare persistence energy to dwell time

rest.DT.PE.r <- readMat(paste(masterdir,'analyses/control_energy/sweep/SweepDTPECorrs_k',numClusters,'.mat',sep = ''))$sweepCorrRest
nback.DT.PE.r <- readMat(paste(masterdir,'analyses/control_energy/sweep/SweepDTPECorrs_k',numClusters,'.mat',sep = ''))$sweepCorrnBack

nperms <- length(rest.DT.PE.r)
grp <- c(rep('Rest',nperms),rep('n-back',nperms))

# compute 1-tailed p-value for hypothesis driven test, where I'm expected rest to be <0 and n-back to >0
pvals <- p.adjust(c(mean(rest.DT.PE.r > 0),mean(nback.DT.PE.r< 0)),method = 'bonf')
p.lab <- pval.label.np(pvals,nperms)
col <- rep('black',2)
col[pvals < 0.05 | pvals > 0.95] <- 'red'

df <- data.frame(r=c(rest.DT.PE.r,nback.DT.PE.r),p=pvals,grp=grp)
p <- ggplot(data=df) + geom_boxplot(aes(y = r,x=grp),color ='#E2492FFF',fill = '#E2492F1A',outlier.size=0.25) +
  scale_color_manual(limits = c('p > 0.05','p < 0.05'),values = c('#71AABE','#005C9F')) + xlab('') + ylab('r(Energy,DT)') +
  scale_y_continuous(limits = c(-1,1)) + scale_x_discrete(limits=c('Rest','n-back')) +
  ggtitle('Min. Persistence Energy') +
  theme_classic() + theme(text = element_text(size = 8)) + theme(plot.title = element_text(size=8,hjust=0.5)) +
  theme(legend.title = element_blank(),legend.key.size = unit(0.5,'line'))
for(Scan in 1:2){
  p <- p + annotate("text", x = c('Rest','n-back')[Scan], y = 1,label = p.lab[Scan],color = col[Scan],size=2)
}
p
ggsave(plot = p,filename = paste(savedir,'PersistEnergyDTCorrSweepBootstrap_k',numClusters,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2,useDingbats=FALSE)