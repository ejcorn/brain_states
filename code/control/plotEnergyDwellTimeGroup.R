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
onDiag <- (1:numClusters) + (numClusters*(0:(numClusters-1)))
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  
clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")

# compare persistence energy to persistence probabilities

energy <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k_',numClusters,'.mat',sep = ''))$Epersist
null_energy <- readMat(paste(masterdir,'analyses/control_energy/DLWNullsPersistenceEnergy_k',numClusters,'.mat',sep = ''))$DLWNullPersistenceEnergy

rdt <- colMeans(readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[1],'DwellTime_k',numClusters,name_root,".mat",sep = ""))$DwellTimeMean)
ndt <- colMeans(readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[2],'DwellTime_k',numClusters,name_root,".mat",sep = ""))$DwellTimeMean)

nperms <- nrow(null_energy)
rest.PP <- lapply(1:nperms, function(P) cor.test(null_energy[P,],rdt))
rest.DT.r <- sapply(rest.PP,function(x) x$estimate)
rest.DT.sig <- as.numeric(sapply(rest.PP,function(x) x$p.value < 0.05))
nback.PP <- lapply(1:nperms, function(P) cor.test(null_energy[P,],ndt))
nback.DT.r <- sapply(nback.PP,function(x) x$estimate)
nback.DT.sig <- as.numeric(sapply(nback.PP,function(x) x$p.value < 0.05))
grp <- c(rep('Rest',nperms),rep('n-back',nperms))
sig <- ifelse(c(rest.DT.sig,nback.DT.sig),yes = 'p < 0.05',no= 'p > 0.05');
actual <- lapply(list(rdt,ndt),function(x) cor.test(x,as.numeric(energy)))
actual <- sapply(actual, function(x) x$estimate)# * (x$p.value < 0.05))
# compute 1-tailed p-value for hypothesis driven test, where I'm expected rest to be lower corr than null and n-back to be higher than null
pvals <- p.adjust(c(mean(actual[1] > rest.DT.r),mean(actual[2] < nback.DT.r)),method = 'bonf')
p.lab <- pval.label.np(pvals,nperms)
col <- rep('black',numClusters)
col[pvals < 0.05 | pvals > 0.95] <- 'red'

p <- ggplot() + geom_boxplot(aes(y = c(rest.DT.r,nback.DT.r),x=grp),color ='#71AABE',fill = '#71AABE1A') +
  geom_point(aes(x=c('Rest','n-back'),y=actual),color = '#CE2B37',fill = '#CE2B37',size = 1,shape = 25) +
  scale_color_manual(limits = c('p > 0.05','p < 0.05'),values = c('#71AABE','#005C9F')) + xlab('') + ylab('r(Energy,DT)') +
  scale_y_continuous(limits = c(-1,1)) + scale_x_discrete(limits=c('Rest','n-back')) +
  ggtitle('Min. Persistence Energy') +
  theme_classic() + theme(text = element_text(size = 8)) + theme(plot.title = element_text(size=8,hjust=0.5)) +
  theme(legend.title = element_blank(),legend.key.size = unit(0.5,'line'))
for(Scan in 1:2){
  p <- p + annotate("text", x = c('Rest','n-back')[Scan], y = 1,label = p.lab[Scan],color = col[Scan],size=2)
}
p
ggsave(plot = p,filename = paste(savedir,'MinEnergyDTCorrVsNull_k',numClusters,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2,useDingbats=FALSE)
