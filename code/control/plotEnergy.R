args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])

name_root <- 'ScanCLaus250Z0corrfinal'
numClusters <- 5
masterdir <- paste("~/Dropbox/Cornblath_Bassett_Projects/code/control_fc/restnbackpipeline/results/",name_root,"/",sep="")
#masterdir <- paste("/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_",name_root,"/",sep="")

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

scanlab <- c("RestComb","nBackComb")
persistExclude = (1:numClusters) + (numClusters*(0:(numClusters-1)))
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
clusterAssignments <- readMat(paste('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterAssignments$clusterAssignments[[1]][[5]])
kClusterCentroids <- unlist(clusterAssignments$clusterAssignments[[1]][[1]])
clusterNames <- c("DMN+", "DMN-", "FPN+", "VIS+", "VIS-")
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")

source('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode/GeomSplitViolin.R')
source('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode/plottingfxns.R')
#source('~/Dropbox/Cornblath_Bassett_Projects/code/control_fc/restnbackpipeline/analysiscode/GeomSplitViolin.R')
#source('~/Dropbox/Cornblath_Bassett_Projects/code/control_fc/restnbackpipeline/analysiscode/plottingfxns.R')

savedir = paste(masterdir,'/analyses/control_energy',sep='');
if(!dir.exists(savedir)){
  dir.create(savedir,recursive = TRUE)
}
# account for distance between centroids

persistExclude = (1:numClusters) + (numClusters*(0:(numClusters-1)))
InterStateDistance = matrix(0,nrow = numClusters,ncol = numClusters)
for(K1 in 1:numClusters){
  for(K2 in 1:numClusters){
    InterStateDistance[K1,K2] = sqrt(sum((kClusterCentroids[,K1] - kClusterCentroids[,K2])^2))
    if(K1 == K2){
      InterStateDistance[K1,K2] = sqrt(sum(kClusterCentroids[,K1]^2))
    }
  }
}

InterStateDistance = as.vector(t(InterStateDistance))
# compare energy to null

energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$E
null_energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$ENull.Rewire
rtp <- colMeans(readMat(paste(masterdir,"analyses/transitionprobabilities/",
	scanlab[1],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability)
ntp <- colMeans(readMat(paste(masterdir,"analyses/transitionprobabilities/",
	scanlab[2],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability)

nperms <- nrow(null_energy)
energy.repmat <- kronecker(matrix(1,nperms,1),energy)
E.gt.null <- colMeans(energy.repmat > null_energy)
E.gt.null <- matrix(E.gt.null,nrow = numClusters,ncol = numClusters,byrow = T)
m <- MatrixPrep(E.gt.null)
cm <- MakeColormap(m)
p <- ggplot() + geom_tile(aes(x=m$Var1,y=m$Var2,fill=m$value)) +
  #scale_fill_gradientn(colours = cm,name ="") +
  scale_fill_viridis(option = 'plasma') + coord_fixed() +
  scale_x_discrete(limits = 1:numClusters, breaks = 1:numClusters, labels = list(clusterNames), expand = c(0,0)) + 
  scale_y_discrete(limits = numClusters:1, breaks = numClusters:1, labels = list(clusterNames), expand = c(0,0)) + 
  theme_bw() + xlab('') + ylab('') + ggtitle('Min. Energy > Null') + theme(plot.title = element_text(size=8,hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(axis.ticks.x = element_blank()) + theme(axis.ticks.y = element_blank()) + 
  theme(text = element_text(size=8)) + theme(legend.position = 'bottom',legend.title = element_blank(),legend.key.size = unit(0.5,'line'))
p

if(numClusters == 5){
  p <- p + theme(axis.text.x = element_text(color = clusterColors),axis.text.y = element_text(color = clusterColors))
}
p

ggsave(plot = p,filename = paste(savedir,'/TransPersistMinEnergyVsNull_k',numClusters,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2)


# compare persistence energy to null

energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$E[persistExclude]
null_energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$ENull.Rewire[,persistExclude]
nperms <- nrow(null_energy)
energy.repmat <- kronecker(matrix(1,nperms,1),matrix(energy,nrow=1))
pvals <- colMeans(energy.repmat > null_energy)
null_energy <- as.vector(null_energy)
states <- as.vector(matrix(clusterNames,nrow=nperms,ncol=numClusters))
col <- rep('black',numClusters)
col[pvals < 0.05 | pvals > 0.95] <- 'red'

p <- ggplot() + geom_boxplot(aes(y = null_energy, x= states),color ='#71AABE',fill = '#71AABE1A') +
  geom_point(aes(x= clusterNames,y=energy),fill = '#CE2B37',color = '#CE2B37',size = 1,shape = 25) +
  geom_text(aes(x = clusterNames, y = 1.1*max(null_energy),label = paste('p =',signif(pvals,2)),color = col),size=2) +
  #scale_color_manual(limits = c('p > 0.05','p < 0.05'),values = c('#71AABE','#005C9F')) + 
  xlab('') + ylab('Min. Control Energy') +
  ggtitle('Persistence Energy') +
  theme_classic() + theme(text = element_text(size = 8)) + theme(plot.title = element_text(size=8,hjust=0.5)) +
  theme(legend.position = 'none')
if(numClusters == 5){
  p <- p + theme(axis.text.x = element_text(color = clusterColors))
}
#p

ggsave(plot = p,filename = paste(savedir,'/MinPersistEnergyVsNull_k',numClusters,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2)

# compare correlation between distance and energy to transition between distance-normalized states for null vs actual
energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$E
null_energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$ENull.Rewire
actual.E.d <- cor(energy[-persistExclude],InterStateDistance[-persistExclude])
null.E.d <- sapply(1:nperms, function(P) cor(null_energy[P,-persistExclude],InterStateDistance[-persistExclude]))
pvals <- mean(actual.E.d > null.E.d)
col <- 'black'
col[pvals < 0.05 | pvals > 0.95] <- 'red'

p <- ggplot() + geom_boxplot(aes(y = null.E.d,x=' '),alpha = 0.6,color ='#71AABE') +
  geom_point(aes(x=' ',y=actual.E.d),fill = '#CE2B37',color = '#CE2B37',size = 1,shape = 25) +
  annotate("text", x = ' ', y = 0,label = paste('p =',signif(pvals,2)),color = col,size=2) +
  #scale_color_manual(limits = c('p > 0.05','p < 0.05'),values = c('#71AABE','#005C9F')) + 
  xlab('') + ylab('r(Energy,Distance)') +
  scale_y_continuous(limits = c(-1,0)) + ggtitle('Min. Transition Energy') +
  theme_classic() + theme(text = element_text(size = 8)) + theme(plot.title = element_text(size=8,hjust=0.5)) +
  theme(legend.title = element_blank(),legend.key.size = unit(0.5,'line'))
p

ggsave(plot = p,filename = paste(savedir,'/MinEnergyDistanceCorrVsNull_k',numClusters,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2)

# compare the correlation between energy and TP for actual A and null, rest and n-back

energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$E
null_energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$ENull.Rewire
rest.TP <- lapply(1:nperms, function(P) cor.test(null_energy[P,-persistExclude],rtp[-persistExclude]))
rest.TP.r <- sapply(rest.TP,function(x) x$estimate)
rest.TP.sig <- as.numeric(sapply(rest.TP,function(x) x$p.value < 0.05))
nback.TP <- lapply(1:nperms, function(P) cor.test(null_energy[P,-persistExclude],ntp[-persistExclude]))
nback.TP.r <- sapply(nback.TP,function(x) x$estimate)
nback.TP.sig <- as.numeric(sapply(nback.TP,function(x) x$p.value < 0.05))
grp <- c(rep('Rest',nperms),rep('n-back',nperms))
sig <- ifelse(c(rest.TP.sig,nback.TP.sig),yes = 'p < 0.05',no= 'p > 0.05');
actual <- lapply(list(rtp[-persistExclude],ntp[-persistExclude]),function(x) cor.test(x,energy[-persistExclude]))
actual <- sapply(actual, function(x) x$estimate)# * (x$p.value < 0.05))

# correlate by rows
idx1 = c(1,5,9,13,17)
idx2 = c(4,8,12,16,20)
for(i in 1:numClusters){ 
  print(cor.test(ntp[-persistExclude][idx1[i]:idx2[i]],energy[-persistExclude][idx1[i]:idx2[i]]))
}

pvals <- c(mean(actual[1] > rest.TP.r),mean(actual[2] > nback.TP.r))
col <- rep('black',numClusters)
col[pvals < 0.05 | pvals > 0.95] <- 'red'

p <- ggplot() + geom_boxplot(aes(y = c(rest.TP.r,nback.TP.r),x=grp),alpha = 0.6,color ='#71AABE') +
  geom_point(aes(x=c('Rest','n-back'),y=actual),fill = '#CE2B37',color = '#CE2B37',size = 1,shape = 25) +
  #scale_color_manual(limits = c('p > 0.05','p < 0.05'),values = c('#71AABE','#005C9F')) + 
  xlab('') + ylab('r(Energy,TP)') +
  scale_x_discrete(limits=c('Rest','n-back')) +
  scale_y_continuous(limits = c(0,1)) + ggtitle('Min. Transition Energy') +
  theme_classic() + theme(text = element_text(size = 8)) + theme(plot.title = element_text(size=8,hjust=0.5)) +
  theme(legend.title = element_blank(),legend.key.size = unit(0.5,'line'))
p
for(Scan in 1:2){
  p <- p + annotate("text", x = c('Rest','n-back')[Scan], y = 1,label = paste('p =',signif(pvals[Scan],2)),color = col[Scan],size=2)
}
p

ggsave(plot = p,filename = paste(savedir,'/MinEnergyTPCorrVsNull_k',numClusters,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2)

# compare correlation between energy and PP for actual A and null, rest and n-back

energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$E
null_energy <- readMat(paste(masterdir,'analyses/control_energy/TransPersistEnergyVsRewire_k_',numClusters,'.mat',sep = ''))$ENull.Rewire
rtp <- colMeans(readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[1],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability)
ntp <- colMeans(readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[2],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability)
nperms <- nrow(null_energy)
rest.PP <- lapply(1:nperms, function(P) cor.test(null_energy[P,persistExclude],rtp[persistExclude]))
rest.PP.r <- sapply(rest.PP,function(x) x$estimate)
rest.PP.sig <- as.numeric(sapply(rest.PP,function(x) x$p.value < 0.05))
nback.PP <- lapply(1:nperms, function(P) cor.test(null_energy[P,persistExclude],ntp[persistExclude]))
nback.PP.r <- sapply(nback.PP,function(x) x$estimate)
nback.PP.sig <- as.numeric(sapply(nback.PP,function(x) x$p.value < 0.05))
grp <- c(rep('Rest',nperms),rep('n-back',nperms))
sig <- ifelse(c(rest.PP.sig,nback.PP.sig),yes = 'p < 0.05',no= 'p > 0.05');
actual <- lapply(list(rtp[persistExclude],ntp[persistExclude]),function(x) cor.test(x,energy[persistExclude]))
actual <- sapply(actual, function(x) x$estimate)# * (x$p.value < 0.05))
pvals <- c(mean(actual[1] > rest.PP.r),mean(actual[2] > nback.PP.r))
col <- rep('black',numClusters)
col[pvals < 0.05 | pvals > 0.95] <- 'red'

p <- ggplot() + geom_boxplot(aes(y = c(rest.PP.r,nback.PP.r),x=grp),color ='#71AABE',fill = '#71AABE1A') +
  geom_point(aes(x=c('Rest','n-back'),y=actual),color = '#CE2B37',fill = '#CE2B37',size = 1,shape = 25) +
  scale_color_manual(limits = c('p > 0.05','p < 0.05'),values = c('#71AABE','#005C9F')) + xlab('') + ylab('r(Energy,PP)') +
  scale_y_continuous(limits = c(-1,1)) + scale_x_discrete(limits=c('Rest','n-back')) +
  ggtitle('Min. Persistence Energy') +
  theme_classic() + theme(text = element_text(size = 8)) + theme(plot.title = element_text(size=8,hjust=0.5)) +
  theme(legend.title = element_blank(),legend.key.size = unit(0.5,'line'))
for(Scan in 1:2){
  p <- p + annotate("text", x = c('Rest','n-back')[Scan], y = 1,label = paste('p =',signif(pvals[Scan],2)),color = col[Scan],size=2)
}
p
ggsave(plot = p,filename = paste(savedir,'/MinEnergyPPCorrVsNull_k',numClusters,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2,useDingbats=FALSE)

