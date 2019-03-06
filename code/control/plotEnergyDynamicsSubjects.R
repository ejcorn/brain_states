# relate subject level control energy metrics to empirically obtained brain state dynamics
# specifically, we expect persistence energy to explain persistence probability and dwell time

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

scanlab <- c("RestComb","nBackComb")
onDiag <- (1:numClusters) + (numClusters*(0:(numClusters-1)))
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  
clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")

# compare persistence energy to persistence probabilities

matResults <- readMat(paste(masterdir,'analyses/control_energy/SubjectPersistenceEnergy_k',numClusters,'.mat',sep = ''))
subjectPersistenceEnergy <- matResults$subjectPersistenceEnergy

rtp <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[1],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability
ntp <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[2],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability

# extract persistence probabilities using indices that correspond to diagonal matrix elements
rpp <- rtp[,onDiag]
npp <- ntp[,onDiag]

nobs <- nrow(rtp)
restPP.cors.within <- sapply(1:nobs, function(N) cor(rpp[N,],as.numeric(subjectPersistenceEnergy[N,])))
nbackPP.cors.within <- sapply(1:nobs, function(N) cor(npp[N,],as.numeric(subjectPersistenceEnergy[N,])))

restPP.cors.across <- lapply(1:numClusters, function(K) cor.test(rpp[,K],subjectPersistenceEnergy[,K],use='pairwise.complete.obs'))
nbackPP.cors.across <- lapply(1:numClusters, function(K) cor.test(npp[,K],subjectPersistenceEnergy[,K],use='pairwise.complete.obs'))

# compare persistence energy to dwell time

rdt <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[1],'DwellTime_k',numClusters,name_root,".mat",sep = ""))$DwellTimeMean
ndt <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[2],'DwellTime_k',numClusters,name_root,".mat",sep = ""))$DwellTimeMean

restDT.cors.within <- sapply(1:nobs, function(N) cor(rdt[N,],as.numeric(subjectPersistenceEnergy[N,])))
nbackDT.cors.within <- sapply(1:nobs, function(N) cor(ndt[N,],as.numeric(subjectPersistenceEnergy[N,])))

restDT.cors.across <- lapply(1:numClusters, function(K) cor.test(rdt[,K],subjectPersistenceEnergy[,K],use='pairwise.complete.obs'))
nbackDT.cors.across <- lapply(1:numClusters, function(K) cor.test(ndt[,K],subjectPersistenceEnergy[,K],use='pairwise.complete.obs'))