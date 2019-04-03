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

# compare persistence energy to switching rate

matResults <- readMat(paste(masterdir,'analyses/control_energy/SubjectPersistenceEnergy_k',numClusters,'.mat',sep = ''))
subjectPersistenceEnergy <- matResults$subjectPersistenceEnergy

# compare mean persistence energy to switch rate

rtn <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[1],'NumTransitions_k',numClusters,name_root,".mat",sep = ""))$numTransitions
ntn <- readMat(paste(masterdir,"analyses/transitionprobabilities/",
                              scanlab[2],'NumTransitions_k',numClusters,name_root,".mat",sep = ""))$numTransitions

# extract persistence probabilities using indices that correspond to diagonal matrix elements
onDiag <- (1:numClusters) + (numClusters*(0:(numClusters-1)))
# if higher persistence energy means stronger propensity to leave state
# then switch rate should be higher for people with higher persistence energy
# what percentage of transitions switch states
r.switchrate <- ( 119 - rowSums(rtn[,onDiag]) )/ 119	# 120 TRs in rest scans, 119 transitions
n.switchrate <- ( 224 - rowSums(ntn[,onDiag]) )/ 224	# 225 TRs in n-back scans, 224 transitions

plot(rowMeans(subjectPersistenceEnergy),r.switchrate)

cor(rowMeans(subjectPersistenceEnergy),n.switchrate)