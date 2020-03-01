args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]
c <- args[4]
nperms <- args[5]

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

# load VIS #

E.VIS <- readMat(paste(masterdir,'analyses/control_energy/model_task_inputs/TransitionEnergies_c',c,'VISWeightedControl_k',
	numClusters,'.mat',sep=''))
T.idx.VIS <- which.min(cor(E.VIS$E.full[offDiag,],nbackTP[offDiag],method=method))
E.T.VIS <- E.VIS$E.full[,T.idx.VIS]
T.opt.VIS <- signif(E.VIS$T.rng[T.idx.VIS],2) # T value at which correlation between 2-back TP and energy is maximal
print(paste('Optimal T for VIS predicting 2-back:',T.opt.VIS))

# now load null models for that T value
E.VIS.DLWNull <- do.call('cbind',lapply(1:nperms, function(P) readMat(paste(masterdir,
	'analyses/control_energy/model_task_inputs/DLWNull_c',c,'T',T.opt.VIS,
	'/DLWNullPerm',P,'_FA_TransitionEnergiesVISWeightedControl_c',c,'T',T.opt.VIS,'_k',numClusters,'.mat',sep=''))$E.full))
E.VIS.BCTNull <- do.call('cbind',lapply(1:nperms, function(P) readMat(paste(masterdir,
	'analyses/control_energy/model_task_inputs/BCTNull_c',c,'T',T.opt.VIS,
	'/BCTNullPerm',P,'_FA_TransitionEnergiesVISWeightedControl_c',c,'T',T.opt.VIS,'_k',numClusters,'.mat',sep=''))$E.full))

### compute p-values ###

# VIS control #

r.RealBrain.rest <- cor(E.T.VIS[offDiag],restTP[offDiag],method=method) # correlation between real brain network energy and rest trans probs
r.DLWNull.rest <- cor(E.VIS.DLWNull[offDiag,],restTP[offDiag],method=method) # correlation between DLW null energy and rest trans probs
p.DLWNull.rest <- pval.label.np(pvals=mean(r.RealBrain.rest >= r.DLWNull.rest),n=length(r.DLWNull.rest))

r.BCTNull.rest <- cor(E.VIS.BCTNull[offDiag,],restTP[offDiag],method=method) # correlation between BCT null energy and rest trans probs
p.BCTNull.rest <- pval.label.np(pvals=mean(r.RealBrain.rest >= r.BCTNull.rest),n=length(r.BCTNull.rest))

r.RealBrain.nback <- cor(E.T.VIS[offDiag],nbackTP[offDiag],method=method) # correlation between real brain network energy and nback trans probs
r.DLWNull.nback <- cor(E.VIS.DLWNull[offDiag,],nbackTP[offDiag],method=method) # correlation between DLW null energy and nback trans probs
p.DLWNull.nback <- pval.label.np(pvals=mean(r.RealBrain.nback >= r.DLWNull.nback),n=length(r.DLWNull.nback))

r.BCTNull.nback <- cor(E.VIS.BCTNull[offDiag,],nbackTP[offDiag],method=method) # correlation between BCT null energy and nback trans probs
p.BCTNull.nback <- pval.label.np(pvals=mean(r.RealBrain.nback >= r.BCTNull.nback),n=length(r.BCTNull.nback))

VIS.pvals.rest <- data.frame(test=c('DLW','BCT'),plabel=c(p.DLWNull.rest,p.BCTNull.rest))
VIS.pvals.nback <- data.frame(test=c('DLW','BCT'),plabel=c(p.DLWNull.nback,p.BCTNull.nback))

### plot energy vs. trans prob ###
# scatter.ETP() is in code/plottingfxns/plottingfxns.R

E.T.VIS <- rank(E.T.VIS[offDiag])
restTP <- rank(restTP[offDiag])
nbackTP <- rank(nbackTP[offDiag])

p.rest <- scatter.ETP(E.T.VIS,restTP,RNcolors[1],ttl='Rest',plabel.df=VIS.pvals.rest,method=method)
save(p.rest,E.T.VIS,restTP,VIS.pvals.rest,file=paste0(savedir,'Fig5d__VISControlRest.RData'))
ggsave(plot = p.rest,filename = paste(savedir,'VISWeightedEnergyVsRestTP_c',c,'T',T.opt.VIS,method,'_k',numClusters,'.pdf',sep =''),
	units = 'in',height = 1.5,width = 1.5)
p.nback <- scatter.ETP(E.T.VIS,nbackTP,RNcolors[2],ttl='2-back',plabel.df=VIS.pvals.nback,method=method)
save(p.rest,E.T.VIS,nbackTP,VIS.pvals.nback,file=paste0(savedir,'Fig5d__VISControl2Back.RData'))
ggsave(plot = p.nback,filename = paste(savedir,'VISWeightedEnergyVs2backTP_c',c,'T',T.opt.VIS,method,'_k',numClusters,'.pdf',sep =''),
	units = 'in',height = 1.5,width = 1.5)
