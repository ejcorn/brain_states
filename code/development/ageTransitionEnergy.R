args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]
c <- args[4]
T <- args[5]

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

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep =""))

#################################
### Age and Transition Energy ###
#################################

savedir <- paste(masterdir,'analyses/development/transprobs/',sep = '')
dir.create(path = savedir,recursive = TRUE)

# whole brain control #

transitionEnergy <- readMat(paste(masterdir,'analyses/control_energy/SubjectTransitionEnergy_FA_c',c,'T',T,'k',
	numClusters,'.mat',sep = ''))$subjectTransitionEnergy

# restTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombTransitionProbabilitiesNoPersist_k',
# 	numClusters,name_root,'.mat',sep = ''))$transitionProbability
# TwoBackBlockTP <- readMat(paste(masterdir,'analyses/nbackblocks/TransProbsNoPersist2back_k',numClusters,name_root,'.mat',sep = ''))$BlockTransitionProbability)

numTransitions <- ncol(transitionEnergy)
transLabels <- as.vector(sapply(1:numClusters, function(i) sapply(1:numClusters,function(j) paste(clusterNames[i],'to',clusterNames[j]))))	#label transitions

TETP.age <- lapply(1:numTransitions, function(T) summary(lm.beta(lm(transitionEnergy[,T] ~ age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
	data = demo))))

TETP.age <- lapply(TETP.age, function(X) X$coefficients[2,c('Standardized','Pr(>|t|)')])

pdat <- cbind(as.data.frame(t(as.data.frame(TETP.age, row.names = c('TEB','TEP')))))
rownames(pdat) <- transLabels
p.list <- list.posthoc.correct(list(pdat$TEP),method = 'bonf')  # bonferroni correct over rest, n-back, 0-2 back
print(p.list)

pdat$TEP <- p.list[[1]] # retrieve corrected p-values for each block/scan

b.mat <- matrix(pdat$TEB, nrow = numClusters, byrow = TRUE)
p.mat <- matrix(pdat$TEP, nrow = numClusters, byrow = TRUE)
p1 <- plot.beta.matrix(b.mat,p.mat,clusterColors,clusterNames,title = 'Whole-Brain Control Energy')
save(p1,b.mat,p.mat,file=paste0(savedir,'Fig6c__TransitionEnergyAge.RData'))
ggsave(plot = p1,filename = paste(savedir,'TransitionEnergyAge_c',c,'T',T,'_k',numClusters,'.pdf',sep =''),units = 'cm',height = 5.5,width = 5.5)