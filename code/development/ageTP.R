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
source(paste(basedir,'code/miscfxns/statfxns.R',sep=''))

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
	data = demo))))
nbackTP.age<- lapply(1:numTransitions, function(T) summary(lm.beta(lm(nbackTP[,T] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
	data = demo))))
names(restTP.age) <- transLabels
names(nbackTP.age) <- transLabels
print(nbackTP.age['DMN- to FPN+'])
print(nbackTP.age['DMN+ to FPN+'])

restTP.age <- lapply(restTP.age, function(X) X$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackTP.age <- lapply(nbackTP.age, function(X) X$coefficients[2,c('Standardized','Pr(>|t|)')])
pdat <- cbind(as.data.frame(t(as.data.frame(restTP.age, row.names = c('RestB','RestP')))),as.data.frame(t(as.data.frame(nbackTP.age, row.names = c('nBackB','nBackP')))))
rownames(pdat) <- transLabels
p.list <- list.posthoc.correct(list(pdat$RestP,pdat$nBackP),method = 'bonf')  # bonferroni correct over rest and n-back
print(p.list)
pdat$RestP <- p.list[[1]] #< 0.05
pdat$nBackP <- p.list[[2]] #< 0.05

b.mat <- matrix(pdat$RestB, nrow = numClusters, byrow = TRUE)
p.mat <- matrix(pdat$RestP, nrow = numClusters, byrow = TRUE)
if(name_root == 'ScanCLaus250Z0final' && numClusters == 5){
	save(b.mat,p.mat, file=paste(savedir,'Fig7e__RestAgeTransProbs_k',numClusters,'.RData',sep =''))
}
p <- plot.beta.matrix(b.mat,p.mat,clusterColors,clusterNames,title = 'Rest')
ggsave(plot = p,filename = paste(savedir,'RestTPAge_k',numClusters,'.pdf',sep =''),units = 'cm',height = 5.5,width = 5.5)

b.mat <- matrix(pdat$nBackB, nrow = numClusters, byrow = TRUE)
p.mat <- matrix(pdat$nBackP, nrow = numClusters, byrow = TRUE)
if(name_root == 'ScanCLaus250Z0final' && numClusters == 5){
	save(b.mat,p.mat, file=paste(savedir,'Fig7f__nBackAgeTransProbs_k',numClusters,'.RData',sep =''))
}
p <- plot.beta.matrix(b.mat,p.mat,clusterColors,clusterNames,title = 'n-back')

ggsave(plot = p,filename = paste(savedir,'nBackTPAge_k',numClusters,'.pdf',sep =''),units = 'cm',height = 5.5,width = 5.5)
