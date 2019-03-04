args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

library(ggplot2)
library(R.matlab)
library(RColorBrewer)
library(lm.beta)
library(reshape2)
library(viridis)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))

clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
RNcolors <- c('#005C9F','#FF8400')  

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep =""))

###########################
### Age and Trans Probs ###
###########################

savedir <- paste(masterdir,'analyses/development/transprobs/',sep = '')
if(!dir.exists(savedir)){
	dir.create(path = savedir,recursive = TRUE)
}

restTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombTransitionProbabilities_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability
nbackTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombTransitionProbabilities_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability

numTransitions <- ncol(restTP)
transLabels <- as.vector(sapply(1:numClusters, function(i) sapply(1:numClusters,function(j) paste(clusterNames[i],'to',clusterNames[j]))))	#label transitions

restTP.age <- lapply(1:numTransitions, function(T) summary(lm.beta(lm(restTP[,T] ~ age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackTP.age<- lapply(1:numTransitions, function(T) summary(lm.beta(lm(nbackTP[,T] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])
names(restTP.age) <- transLabels
names(nbackTP.age) <- transLabels
pdat <- as.data.frame(t(as.data.frame(restTP.age, row.names = c('RestB','RestP'))))
pdat$RestP <- p.adjust(pdat$RestP,method = 'bonf') < 0.05
v <- matrix((pdat$RestB*pdat$RestP), nrow = numClusters, byrow = TRUE)
p <- TP.beta.plot(v,clusterColors,title = 'Rest: Age')
ggsave(plot = p,filename = paste(savedir,'RestTPAge_k',numClusters,'.pdf',sep =''),units = 'cm',height = 6,width = 6)

pdat <- as.data.frame(t(as.data.frame(nbackTP.age, row.names = c('nBackB','nBackP'))))
pdat$nBackP <- p.adjust(pdat$nBackP,method = 'bonf') < 0.05
v <- matrix((pdat$nBackB*pdat$nBackP), nrow = numClusters, byrow = TRUE)
p <- TP.beta.plot(v,clusterColors,title = 'n-back: Age')
ggsave(plot = p,filename = paste(savedir,'nBackTPAge_k',numClusters,'.pdf',sep =''),units = 'cm',height = 6,width = 6)


# restTP.acc <- lapply(1:numTransitions, function(T) summary(lm.beta(lm(Overall_Accuracy ~ restTP[,T] + age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
# 	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])
# nbackTP.acc <- lapply(1:numTransitions, function(T) summary(lm.beta(lm(Overall_Accuracy ~ nbackTP[,T] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
# 	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# pdat <- as.data.frame(t(as.data.frame(restTP.acc, row.names = c('RestB','RestP'))))
# rownames(pdat) <- transLabels
# pdat$RestP <- p.adjust(pdat$RestP,method = 'bonf') < 0.05
# v <- matrix((pdat$RestB*pdat$RestP), nrow = numClusters, byrow = TRUE)
# p <- TP.beta.plot(v,clusterColors,title = 'Rest: Overall Accuracy')
# ggsave(plot = p,filename = paste(savedir,'RestTPOverallAccuracy_k',numClusters,'.pdf',sep =''),units = 'in',height = 2,width = 2)

# pdat <- as.data.frame(t(as.data.frame(nbackTP.acc, row.names = c('nBackB','nBackP'))))
# rownames(pdat) <- transLabels
# pdat$nBackP <- p.adjust(pdat$nBackP,method = 'bonf') < 0.05
# v <- matrix((pdat$nBackB*pdat$nBackP), nrow = numClusters, byrow = TRUE)
# p <- TP.beta.plot(v,clusterColors,title = 'n-back: Overall Accuracy')
# ggsave(plot = p,filename = paste(savedir,'nBackTPOverallAccuracy_k',numClusters,'.pdf',sep =''),units = 'in',height = 2,width = 2)


