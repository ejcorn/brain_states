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

# compare transition energy to transition probabilities
# should be anticorrelated

matResults <- readMat(paste(masterdir,'analyses/control_energy/SubjectTransitionEnergy_k',numClusters,'.mat',sep = ''))
subjectTransitionEnergy <- matResults$subjectTransitionEnergy

rtp <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[1],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability
ntp <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[2],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability

numTrans <- numClusters^2
rtp.e.test <- lapply(1:numTrans, function(T) cor.test(rtp[,T],subjectTransitionEnergy[,T]))
rtp.e.r <- sapply(1:numTrans, function(T) rtp.e.test[[T]]$estimate)
rtp.e.p <- p.adjust(sapply(1:numTrans, function(T) rtp.e.test[[T]]$p.value))

ntp.e.test <- lapply(1:numTrans, function(T) cor.test(ntp[,T],subjectTransitionEnergy[,T]))
ntp.e.r <- sapply(1:numTrans, function(T) ntp.e.test[[T]]$estimate)
ntp.e.p <- p.adjust(sapply(1:numTrans, function(T) ntp.e.test[[T]]$p.value))

# compare sum of transition energies into each state with fractional occupancy

rfo <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[1],'FractionalOccupancy_k',numClusters,name_root,".mat",sep = ""))$FractionalOccupancy
nfo <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[2],'FractionalOccupancy_k',numClusters,name_root,".mat",sep = ""))$FractionalOccupancy

# each transition labeled sequentially
idx.mat <- t(matrix(1:numTrans,nrow=numClusters,ncol=numClusters))
# get energy of transitions ending in each state, including persistence
idx.trans.in <- lapply(1:numClusters, function(K) as.numeric(idx.mat[,K]))

# predicted fractional occupancy based on transition energies
pred.fo <- sapply(1:numClusters, function(K) rowSums(subjectTransitionEnergy[,idx.trans.in[[K]]]))

# should be negative, i.e. more total energy needed to get into state = less occupancy
# if at rest you just follow energy
fo.energy.rest <- lapply(1:numClusters, function(K) cor(pred.fo[,K],rfo[,K]))
# for n-back should be opposite or orthogonal
fo.energy.nback <- lapply(1:numClusters, function(K) cor(pred.fo[,K],nfo[,K]))
