args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]
c <- args[4]

library(ggplot2)
library(R.matlab)
library(RColorBrewer)
library(lm.beta)
library(reshape2)
library(viridis)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))
source(paste(basedir,'code/miscfxns/statfxns.R',sep=''))

clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- getClusterColors(numClusters)
RNcolors <- c('#005C9F','#FF8400')  

savedir <- paste(masterdir,'analyses/control_energy/EvsTP_plots/',sep = '')
dir.create(path = savedir,recursive = TRUE)

method = 'spearman'

### load transition probabilities ###

restTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombTransitionProbabilitiesNoPersist_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability
restTP <- colMeans(restTP) # group average
nbackTP <- readMat(paste(masterdir,'analyses/nbackblocks/TransProbsNoPersist2back_k',numClusters,name_root,'.mat',sep = ''))$BlockTransitionProbability
nbackTP <- colMeans(nbackTP,na.rm=T)

onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1))) # construct indices to select only persistence probabilites
offDiag = 1:(numClusters^2); offDiag <- offDiag[-onDiag] # construct indices to select only transition probabilites

### load control energies ###

E.wholebrain <- readMat(paste(masterdir,'analyses/control_energy/Tsweep/c',c,'_k',numClusters,
    '/GroupAverageTransitionEnergies_k',numClusters,'.mat',sep=''))
#T.idx <- which.min(abs(E.wholebrain$T.rng-T)) # get index of T value to use based on input
T.idx <- which.min(cor(E.wholebrain$E.full[offDiag,],restTP[offDiag],method=method))
T.opt.WB <- signif(E.wholebrain$T.rng[T.idx],2)
E.T.WB <- E.wholebrain$E.full[,T.idx]
E.wholebrain.BCTNull <- E.wholebrain$E.BCTnull[,,T.idx]

# this came out to 5.5 .. even though 5.5 gives max correlation at rest
# can stil use 5 for null testing, conservatively as it's not optimal
T.opt.WB <- 5


E.wholebrain.DLWNull <- readMat(paste(masterdir,'analyses/control_energy/DLWNullsTransitionEnergy_FA_k',
	numClusters,'c',c,'T',T.opt.WB,'.mat',sep=''))$DLWNullTransitionEnergy

### compute p-values ###

# whole brain control #
r.RealBrain.rest <- cor(E.T.WB[offDiag],restTP[offDiag],method=method) # correlation between real brain network energy and rest trans probs
r.DLWNull.rest <- cor(t(E.wholebrain.DLWNull[,offDiag]),restTP[offDiag],method=method) # correlation between DLW null energy and rest trans probs
p.DLWNull.rest <- pval.label.np(pvals=mean(r.RealBrain.rest >= r.DLWNull.rest),n=length(r.DLWNull.rest))

r.BCTNull.rest <- cor(E.wholebrain.BCTNull[offDiag,],restTP[offDiag],method=method) # correlation between BCT null energy and rest trans probs
p.BCTNull.rest <- pval.label.np(pvals=mean(r.RealBrain.rest >= r.BCTNull.rest),n=length(r.BCTNull.rest))

r.RealBrain.nback <- cor(E.T.WB[offDiag],nbackTP[offDiag],method=method) # correlation between real brain network energy and nback trans probs
r.DLWNull.nback <- cor(t(E.wholebrain.DLWNull[,offDiag]),nbackTP[offDiag],method=method) # correlation between DLW null energy and nback trans probs
p.DLWNull.nback <- pval.label.np(pvals=mean(r.RealBrain.nback >= r.DLWNull.nback),n=length(r.DLWNull.nback))

r.BCTNull.nback <- cor(E.wholebrain.BCTNull[offDiag,],nbackTP[offDiag],method=method) # correlation between BCT null energy and nback trans probs
p.BCTNull.nback <- pval.label.np(pvals=mean(r.RealBrain.nback >= r.BCTNull.nback),n=length(r.BCTNull.nback))

wholebrain.pvals.rest <- data.frame(test=c('DLW','BCT'),plabel=c(p.DLWNull.rest,p.BCTNull.rest))
wholebrain.pvals.nback <- data.frame(test=c('DLW','BCT'),plabel=c(p.DLWNull.nback,p.BCTNull.nback))

### plot energy vs. trans prob ###
# scatter.ETP() is in code/plottingfxns/plottingfxns.R

# convert to rank values so linear fit makes sense 
E.T.WB <- rank(E.T.WB[offDiag])
restTP <- rank(restTP[offDiag])
nbackTP <- rank(nbackTP[offDiag])

init <- rep(1:numClusters,each=numClusters-1) # starting cluster for each off-diagonal transition
term <- rep(1:numClusters,numClusters)[offDiag] # ending cluster for each off-diagonal transition

# note: rest looks like there are actualy several upward trends that on average form a downward trend
# I checked by coloring the points with init and term above, which revealed that
# these lines do not just reflect rows or columns of the transition matrix
p.rest <- scatter.ETP(E.T.WB,restTP,RNcolors[1],ttl='Rest',plabel.df=wholebrain.pvals.rest,method=method)
save(p.rest,E.T.WB,restTP,wholebrain.pvals.rest,file=paste0(savedir,'Fig5c__WholeBrainControlRest.RData'))
ggsave(plot = p.rest,filename = paste(savedir,'WholeBrainEnergyVsRestTP_c',c,'T',T.opt.WB,method,'_k',numClusters,'.pdf',sep =''),
	units = 'in',height = 1.5,width = 1.5)
p.nback <- scatter.ETP(E.T.WB,nbackTP,RNcolors[2],ttl='2-back',plabel.df=wholebrain.pvals.nback,method=method)
save(p.nback,E.T.WB,nbackTP,wholebrain.pvals.nback,file=paste0(savedir,'Fig5c__WholeBrainControl2Back.RData'))
ggsave(plot = p.nback,filename = paste(savedir,'WholeBrainEnergyVs2backTP_c',c,'T',T.opt.WB,method,'_k',numClusters,'.pdf',sep =''),
	units = 'in',height = 1.5,width = 1.5)
