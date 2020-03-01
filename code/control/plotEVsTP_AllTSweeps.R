# plot weighted control energy vs. transition probabilities for a range of T values
# compare to correlation between system-weighted interstate distance and transition probabilities
# do this for the 7 Yeo cognitive systems
# then also for whole brain control, with comparison to null

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
library(gridExtra)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))
source(paste(basedir,'code/miscfxns/statfxns.R',sep=''))

clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- getClusterColors(numClusters)
RNcolors <- c('#005C9F','#FF8400')  

savedir <- paste(masterdir,'analyses/control_energy/EvsTP_plots/WeightedControl/',sep = '')
dir.create(path = savedir,recursive = TRUE)

### load transition probabilities ###

restTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombTransitionProbabilitiesNoPersist_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability
restTP <- colMeans(restTP) # group average
nbackTP <- readMat(paste(masterdir,'analyses/nbackblocks/TransProbsNoPersist2back_k',numClusters,name_root,'.mat',sep = ''))$BlockTransitionProbability
nbackTP <- colMeans(nbackTP,na.rm=T)

onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1))) # construct indices to select only persistence probabilites
offDiag = 1:(numClusters^2); offDiag <- offDiag[-onDiag] # construct indices to select only transition probabilites

Centroids <- readMat(paste(masterdir,'analyses/centroids/OverallClusterCentroids_k',
	numClusters,name_root,'.mat',sep = ''))$kClusterCentroids
xf <- Centroids[,rep(1:numClusters,times=numClusters)]
x0 <- Centroids[,rep(1:numClusters,each=numClusters)]
InterStateDistance <- colSums((xf - x0)^2)

### load weighted control energies ###

yeo <- readMat(paste(basedir,'data/yeo7netlabelsLaus250.mat',sep=''))$network7labels
yeo <- yeo[-length(yeo)] # remove brainstem
Systems <- c('VIS','SOM','DAT','VAT','LIM','FPN','DMN')
B.matrix <- lapply(1:length(Systems), function(S) sapply(1:numClusters^2, function(K) (1+(yeo==S))^2)) # for weighted control, reconstruct B matrix to use for weighted distance calculation
names(B.matrix) <- Systems
InterStateDistance.weighted <- lapply(Systems, function(S) colSums((xf-x0)^2/B.matrix[[S]])) # interstate distance weighted by input matrix
names(InterStateDistance.weighted) <- Systems

weightedcontrol <- lapply(Systems, function(Sys) readMat(paste(masterdir,'analyses/control_energy/model_task_inputs/TransitionEnergies_c',c,Sys,'WeightedControl_k',
	numClusters,'.mat',sep='')))
names(weightedcontrol) <- Systems
E.null <- weightedcontrol[['SOM']]$E.identity[offDiag,] # get whole brain control from a random system, doesn't matter bc E.identity is always the same. 
E.null <- array(E.null,dim=c(dim(E.null)[1],1,dim(E.null)[2])) # add a 3rd dimension
E.null = list(WB=E.null)

### load whole brain control energies ###

E.wholebrain <- readMat(paste(masterdir,'analyses/control_energy/Tsweep/c',c,'_k',numClusters,
    '/GroupAverageTransitionEnergies_k',numClusters,'.mat',sep=''))

E.brain <- E.wholebrain$E.full[offDiag,]
E.null.full <- list(BCT=E.wholebrain$E.BCTnull[offDiag,,])
t.rng <- as.numeric(E.wholebrain$T.rng)

for(method in c('pearson','spearman')){
	rest.plots.weighted <- lapply(Systems, function(S) p.Tsweep(as.numeric(weightedcontrol[[S]]$T.rng),weightedcontrol[[S]]$E.full[offDiag,],E.null,restTP[offDiag],
		InterStateDistance.weighted[[S]][offDiag],RNcolors[1],ttl=paste('Rest: ',S,'-Weighted',sep=''),method=method,leg=FALSE))
	# tack on whole brain control plot for rest
	rest.plots.uniform <- list(p.Tsweep(t.rng,E.brain,E.null.full,restTP[offDiag],InterStateDistance[offDiag],RNcolors[1],ttl='Rest: Uniform',method=method,leg = FALSE))
	rest.plots <- c(rest.plots.uniform,rest.plots.weighted)
	ggsave(plot = arrangeGrob(grobs=rest.plots,nrow=2),filename = paste(savedir,'AllTSweepsControlEnergyVsRestTP_c',c,method,'_k',numClusters,'.pdf',sep =''),
		units = 'cm',height = 9,width = 18)
	nback.plots.weighted <- lapply(Systems, function(S) p.Tsweep(as.numeric(weightedcontrol[[S]]$T.rng),weightedcontrol[[S]]$E.full[offDiag,],E.null,nbackTP[offDiag],
		InterStateDistance.weighted[[S]][offDiag],RNcolors[2],ttl=paste('2-back: ',S,'-Weighted',sep=''),method=method,leg=FALSE))
	# tack on whole brain control plot for nback
	nback.plots.uniform <- list(p.Tsweep(t.rng,E.brain,E.null.full,nbackTP[offDiag],InterStateDistance[offDiag],RNcolors[2],ttl='2-back: Uniform',method=method,leg = FALSE))
	nback.plots <- c(nback.plots.uniform,nback.plots.weighted)
	ggsave(plot = arrangeGrob(grobs=nback.plots,nrow=2),filename = paste(savedir,'AllTSweepsControlEnergyVs2BackTP_c',c,method,'_k',numClusters,'.pdf',sep =''),
		units = 'cm',height = 9,width = 18)
}

