# relate group level control energy metrics to empirically obtained brain state dynamics
# namely, we expect persistence energy to explain persistence probability and dwell time

args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]
c <- as.numeric(args[4])
T <- as.numeric(args[5])

# variables that get passed in:
# c: normalization factor
# normalize each matrix by 1/c+max(eig(A)), then subtract identity. ensures max eigenvalue is [0,-1), such that matrices are stable and 
# at long time horizons all nodes either go to 0 or to a constant non-zero eigenvector associated with max eig
# T: control horizon, time over which to exert control

masterdir <- paste(basedir,'results/',name_root,'/',sep='')
savedir <- paste(masterdir,'analyses/control_energy/',sep='')

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))
source(paste(basedir,'code/miscfxns/statfxns.R',sep=''))

scanlab <- c("RestComb","nBackComb")
persistExclude = (1:numClusters) + (numClusters*(0:(numClusters-1)))
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  
clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")

# compare persistence energy to null, 

energy <- as.numeric(readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k',numClusters,'c',c,'T',T,'.mat',sep = ''))$Epersist)
null_energy <- readMat(paste(masterdir,'analyses/control_energy/DLWNullsPersistenceEnergy_k',numClusters,'c',c,'T',T,'.mat',sep = ''))$DLWNullPersistenceEnergy

nperms <- nrow(null_energy)
pvals <- p.adjust(sapply(1:numClusters, function(k) mean(energy[k] > null_energy[,k])),method='bonferroni')		# compute 1-tailed p-values
null_energy <- as.vector(null_energy)
states <- as.vector(matrix(clusterNames,nrow=nperms,ncol=numClusters))
col <- rep('black',numClusters)
col[pvals < 0.05 | pvals > 0.95] <- 'red'
p.lab <- pval.label.np(pvals,nperms)		# make p-value labels that don't say p = 0

p <- ggplot() + geom_boxplot(aes(y = null_energy, x= states),color ='#71AABE',fill = '#71AABE1A') +
  geom_point(aes(x= clusterNames,y=energy),fill = '#CE2B37',color = '#CE2B37',size = 1,shape = 25) +
  geom_text(aes(x = clusterNames, y = 1.1*max(null_energy),label = p.lab,color=col),size=1.5) +
  scale_x_discrete(limits = clusterNames, breaks=clusterNames, labels = list(clusterNames)) +
  scale_color_manual(limits = c('black','red'),values = c('black','red')) + 
  xlab('') + ylab('Min. Control Energy') +
  ggtitle('Persistence Energy') +
  theme_classic() + theme(text = element_text(size = 8)) + theme(plot.title = element_text(size=8,hjust=0.5)) +
  theme(legend.position = 'none')
if(numClusters == 5){
  p <- p + theme(axis.text.x = element_text(color = clusterColors))
}
#p

ggsave(plot = p,filename = paste(savedir,'MinPersistEnergyVsNull_k',numClusters,'c',c,'T',T,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2)