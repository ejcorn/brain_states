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
clusterColors <- getClusterColors(numClusters)
RNcolors <- c('#005C9F','#FF8400')  

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep =""))

###########################
### Age and Trans Probs ###
###########################

savedir <- paste(masterdir,'analyses/development/transprobs/',sep = '')
dir.create(path = savedir,recursive = TRUE)

restTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombTransitionProbabilitiesNoPersist_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability
nbackTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombTransitionProbabilitiesNoPersist_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability
nbackBlockTP <- lapply(0:2, function(B) readMat(paste(masterdir,'analyses/nbackblocks/TransProbsNoPersist',B,'back_k',numClusters,name_root,'.mat',sep = ''))$BlockTransitionProbability)

numTransitions <- ncol(restTP)
transLabels <- as.vector(sapply(1:numClusters, function(i) sapply(1:numClusters,function(j) paste(clusterNames[i],'to',clusterNames[j]))))	#label transitions

restTP.age <- lapply(1:numTransitions, function(T) summary(lm.beta(lm(restTP[,T] ~ age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
	data = demo))))
nbackTP.age<- lapply(1:numTransitions, function(T) summary(lm.beta(lm(nbackTP[,T] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
	data = demo))))
nbackBlockTP.age<- lapply(nbackBlockTP, function(blockTP) lapply(1:numTransitions, function(T) summary(lm.beta(lm(blockTP[,T] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
	data = demo)))))
names(restTP.age) <- transLabels
names(nbackTP.age) <- transLabels

for(B in 1:3){names(nbackBlockTP.age[[B]]) <- transLabels}

restTP.age <- lapply(restTP.age, function(X) X$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackTP.age <- lapply(nbackTP.age, function(X) X$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackBlockTP.age <- lapply(nbackBlockTP.age, function(X) lapply(X, function(X) X$coefficients[2,c('Standardized','Pr(>|t|)')]))

pdat <- cbind(as.data.frame(t(as.data.frame(restTP.age, row.names = c('RestB','RestP')))),as.data.frame(t(as.data.frame(nbackTP.age, row.names = c('nBackB','nBackP')))))
rownames(pdat) <- transLabels
pdat.blocks <- do.call('cbind',lapply(1:3, function(B) as.data.frame(t(as.data.frame(nbackBlockTP.age[[B]], row.names = c(paste(B-1,'BackB',sep=''),paste(B-1,'BackP',sep='')))))))
rownames(pdat.blocks) <- transLabels
pdat.all <- cbind(pdat,pdat.blocks)
p.list <- list.posthoc.correct(list(pdat.all$RestP,pdat.all$nBackP,pdat.all$`0BackP`,pdat.all$`1BackP`,pdat.all$`2BackP`),method = 'bonf')  # bonferroni correct over rest, n-back, 0-2 back
print(p.list)

pdat.all$RestP <- p.list[[1]] # retrieve corrected p-values for each block/scan
pdat.all$nBackP <- p.list[[2]]
pdat.all$`0BackP` <- p.list[[3]]
pdat.all$`1BackP` <- p.list[[4]]
pdat.all$`2BackP` <- p.list[[5]]

b.list <- list(pdat.all$RestB,pdat.all$nBackB,lapply(0:2, function(B) get(paste(B,'BackB',sep=''),pdat.all))) # get all betas in one list
bmax <- max(unlist(b.list),na.rm=T) # get max beta value across all matrices
bmin <- min(unlist(b.list),na.rm=T) # get min beta value across all matrices

b.mat <- matrix(pdat.all$RestB, nrow = numClusters, byrow = TRUE)
p.mat <- matrix(pdat.all$RestP, nrow = numClusters, byrow = TRUE)
p1 <- plot.beta.matrix(b.mat,p.mat,clusterColors,clusterNames,title = 'Rest',bmin,bmax)
save(p1,b.mat,p.mat,file=paste0(savedir,'Fig6d__RestTPAge.RData'))
ggsave(plot = p1,filename = paste(savedir,'RestTPAge_k',numClusters,'.pdf',sep =''),units = 'cm',height = 5.5,width = 5.5)

b.mat <- matrix(pdat.all$nBackB, nrow = numClusters, byrow = TRUE)
p.mat <- matrix(pdat.all$nBackP, nrow = numClusters, byrow = TRUE)
p2 <- plot.beta.matrix(b.mat,p.mat,clusterColors,clusterNames,title = 'n-back',bmin,bmax)
save(p2,b.mat,p.mat,file=paste0(savedir,'Fig6h__nbackTPAge.RData'))
ggsave(plot = p2,filename = paste(savedir,'nBackTPAge_k',numClusters,'.pdf',sep =''),units = 'cm',height = 5.5,width = 5.5)

b.mat <- matrix(pdat.all$`0BackB`, nrow = numClusters, byrow = TRUE)
p.mat <- matrix(pdat.all$`0BackP`, nrow = numClusters, byrow = TRUE)
p3 <- plot.beta.matrix(b.mat,p.mat,clusterColors,clusterNames,title = '0-back',bmin,bmax)
save(p3,b.mat,p.mat,file=paste0(savedir,'Fig6e__0backTPAge.RData'))
ggsave(plot = p3,filename = paste(savedir,'0BackTPAge_k',numClusters,'.pdf',sep =''),units = 'cm',height = 5.5,width = 5.5)

b.mat <- matrix(pdat.all$`1BackB`, nrow = numClusters, byrow = TRUE)
p.mat <- matrix(pdat.all$`1BackP`, nrow = numClusters, byrow = TRUE)
p4 <- plot.beta.matrix(b.mat,p.mat,clusterColors,clusterNames,title = '1-back',bmin,bmax)
save(p4,b.mat,p.mat,file=paste0(savedir,'Fig6f__1backTPAge.RData'))
ggsave(plot = p4,filename = paste(savedir,'1BackTPAge_k',numClusters,'.pdf',sep =''),units = 'cm',height = 5.5,width = 5.5)

b.mat <- matrix(pdat.all$`2BackB`, nrow = numClusters, byrow = TRUE)
p.mat <- matrix(pdat.all$`2BackP`, nrow = numClusters, byrow = TRUE)
p5 <- plot.beta.matrix(b.mat,p.mat,clusterColors,clusterNames,title = '2-back',bmin,bmax)
save(p5,b.mat,p.mat,file=paste0(savedir,'Fig6g__2backTPAge.RData'))
ggsave(plot = p5,filename = paste(savedir,'2BackTPAge_k',numClusters,'.pdf',sep =''),units = 'cm',height = 5.5,width = 5.5)

