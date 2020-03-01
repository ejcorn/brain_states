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

nbbeh <- read.csv(paste(basedir,'data/n1601_nbackBehavior_from_20160207_dataRelease.csv',sep=''))
nbbeh <- nbbeh[which(nbbeh$scanid %in% demo$scanid),]
nbbeh <- nbbeh[order(nbbeh$scanid),]

cnb <- read.csv(paste(basedir,'data/n1601_cnb_factor_scores_tymoore_20151006.csv',sep=''))
cnb <- cnb[which(nbbeh$scanid %in% demo$scanid),]
cnb <- cnb[order(cnb$scanid),]

data <- cbind(demo, cnb[,grepl('Accuracy',colnames(cnb))], nbbeh[,grepl('Dprime',colnames(nbbeh))])

# can control for GLM activation betas and transitions still explain behavior better
# betas <- read.csv(paste(basedir,'data/ELIDANI_lausannescale125.csv',sep=''))
# betas <- betas[which(betas$scanid %in% demo$scanid),]
# betas <- betas[order(betas$scanid),]
# yeo <- readMat(paste(basedir,'data/yeo7netlabelsLaus125.mat',sep=''))$network7labels
# # select fpn betas only
# betas.fpn <- betas[,paste("Laussannescale125_contrast4_2back.0back_roi",which(yeo==6),sep='')]
# #betas.fpn <- data.frame(meanfpnbeta = rowMeans(betas.fpn))
# betas.fpn <- cbind(data.frame(scanid=betas$scanid),betas.fpn)

# data <- merge(data,betas.fpn,by='scanid',all=TRUE)

#################################
### Cognition and Trans Probs ###
#################################

savedir <- paste(masterdir,'analyses/cognition/transprobs/',sep = '')
dir.create(path = savedir,recursive = TRUE)

nbackTP <- readMat(paste(masterdir,'analyses/nbackblocks/TransProbsNoPersist2back_k',numClusters,name_root,'.mat',sep = ''))$BlockTransitionProbability

numTransitions <- ncol(nbackTP)
transLabels <- as.vector(sapply(1:numClusters, function(i) sapply(1:numClusters,function(j) paste(clusterNames[i],'to',clusterNames[j]))))	#label transitions
colnames(nbackTP) <- transLabels

dprime <- data$nbackBehTwobackDprime
data <- cbind(data,dprime,nbackTP) # concatenate data

# make a new data frame for each transition
covariates <- c('dprime','age_in_yrs', 'BrainSegVol' ,'handedness' ,'nbackRelMeanRMSMotion' ,'Sex')

nbackTP.dprime.data <- lapply(transLabels, function(T) data[,c(T,covariates)])
names(nbackTP.dprime.data) <- transLabels

# regress dprime on all remaining variables in transition-specific dataframe, extract coefficients and p-values
nbackTP.dprime.coefs <- list()
nbackTP.dprime <- list()

for(T in transLabels){
	nbackTP.dprime.coefs[[T]] <- summary(lm.beta(lm(dprime ~ ., 
	data = nbackTP.dprime.data[[T]])))$coefficients
	
	# if it's a persistence prob, save NAs for beta and p-value because they're all 0
	if(length(nbackTP.dprime.coefs[[T]]) == (5*(length(covariates)+1))){
		nbackTP.dprime[[T]] <- nbackTP.dprime.coefs[[T]][paste('`',T,'`',sep=''),c('Standardized','Pr(>|t|)')]
	} else if(length(nbackTP.dprime.coefs[[T]]) < (5*(length(covariates)+1))){
		# R will automatically delete the row of coefficients for transition prob because it's singular
		# I can detect this by checking the size of the output coefficient matrix
		nbackTP.dprime[[T]] <- c(NA,NA)
	}
	
}

pdat <- cbind(as.data.frame(as.data.frame(t(as.data.frame(nbackTP.dprime, row.names = c('nBackB','nBackP'))))))
rownames(pdat) <- transLabels
p.list <- list.posthoc.correct(list(pdat$nBackP),method = 'bonf')  # bonferroni correct over rest and n-back
pdat$nBackP <- p.list[[1]]

# test reshaping to matrix:
# matrix(1:numClusters, nrow = numClusters, byrow = TRUE)

b.mat <- matrix(pdat$nBackB, nrow = numClusters, byrow = TRUE)
p.mat <- matrix(pdat$nBackP, nrow = numClusters, byrow = TRUE)
p <- plot.beta.matrix(b.mat,p.mat,clusterColors,clusterNames,title = '2-back: d-prime')
ggsave(plot = p,filename = paste(savedir,'TwoBackTPNoPersistDprime_k',numClusters,'.pdf',sep =''),units = 'cm',height = 6,width = 6)
writeMat(paste(savedir,'TwoBackTPNoPersistDprime_k',numClusters,'.mat',sep =''),bmat=b.mat,pmat=p.mat)