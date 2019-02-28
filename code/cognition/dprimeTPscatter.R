args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])

library(ggplot2)
library(R.matlab)
library(RColorBrewer)
library(lm.beta)
library(reshape2)
library(viridis)

masterdir <- paste("/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_",name_root,"/",sep="")
source('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode/plottingfxns.R')

clusterNames <- readMat(paste('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
RNcolors <- c('#005C9F','#FF8400') 

demo <- read.csv(paste('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/Demographics',name_root,'.csv',sep =""))

nbbeh <- read.csv('/data/jag/bassett-lab/Lausanne1601/nback/n1601_nbackBehavior_from_20160207_dataRelease.csv')
nbbeh <- nbbeh[which(nbbeh$scanid %in% demo$scanid),]
nbbeh <- nbbeh[order(nbbeh$scanid),]

cnb <- read.csv('/data/jag/bassett-lab/Lausanne1601/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv')
cnb <- cnb[which(nbbeh$scanid %in% demo$scanid),]
cnb <- cnb[order(cnb$scanid),]

data <- cbind(demo, cnb[,grepl('Accuracy',colnames(cnb))], nbbeh[,grepl('Dprime',colnames(nbbeh))])

####################################
### plot significant transitions ###
####################################

savedir <- paste(masterdir,'analyses/cognition/transprobs/',sep = '')
if(!dir.exists(savedir)){
	dir.create(path = savedir)
}

restTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombTransitionProbabilities_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability
nbackTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombTransitionProbabilities_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability

numTransitions <- ncol(restTP)
transLabels <- as.vector(sapply(1:numClusters, function(i) sapply(1:numClusters,function(j) paste(clusterNames[i],'to',clusterNames[j]))))	#label transitions

# overall d prime
restTP.dprime <- lapply(1:numTransitions, function(T) summary(lm.beta(lm(nbackBehAllDprime ~ restTP[,T] + age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
	data = data)))$coefficients[2,c('Estimate','Pr(>|t|)')])

pdat <- as.data.frame(t(as.data.frame(restTP.dprime, row.names = c('RestB','RestP'))))
rownames(pdat) <- 1:numTransitions
pdat$RestPCorrect <- p.adjust(pdat$RestP,method = 'bonf') < 0.05
transInd <- which(pdat$RestPCorrect)[1]
dp.resid <- residuals(lm(nbackBehAllDprime ~ restTP[,transInd] + age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex,data=data))
d <- as.data.frame(cbind(x = restTP[,transInd],y= dp.resid + restTP[,transInd]*pdat$RestB[transInd]))

p <- ggplot(data = d,aes(x=x,y=y)) + geom_point(color = 'grey70',stroke =0,size = 1, alpha = 0.6) + 
	geom_smooth(fill=RNcolors[1],color=RNcolors[1],method = 'lm') + xlab(paste('P(',transLabels[transInd],')',sep='')) +
	geom_text(aes(x=0.06+ min(d$x),y=max(d$y),label = paste('p = ',signif(pdat$RestP[transInd],2),sep='')),size=2.5) +
	ggtitle('Rest') + ylab('WM Performance') + theme_classic() + theme(text = element_text(size = 8)) +
	theme(plot.title = element_text(size = 8,hjust = 0.5))
ggsave(plot = p,filename = paste(savedir,'Rest',transLabels[transInd],'Dprime_k',numClusters,'.pdf',sep =''),units = 'in',height = 2,width = 2.5)

nbackTP.dprime <- lapply(1:numTransitions, function(T) summary(lm.beta(lm(nbackBehAllDprime ~ nbackTP[,T] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
	data = data)))$coefficients[2,c('Estimate','Pr(>|t|)')])

pdat <- as.data.frame(t(as.data.frame(nbackTP.dprime, row.names = c('nBackB','nBackP'))))
rownames(pdat) <- 1:numTransitions
pdat$nBackPCorrect <- p.adjust(pdat$nBackP,method = 'bonf') < 0.05
transInd <- which(pdat$nBackPCorrect)[1]
dp.resid <- residuals(lm(nbackBehAllDprime ~ nbackTP[,transInd] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,data=data))
d <- as.data.frame(cbind(x = nbackTP[,transInd], y = dp.resid + nbackTP[,transInd]*pdat$nBackB[transInd]))		# store partial residuals for plotting

p <- ggplot(data = d,aes(x=x,y=y)) + geom_point(color = 'grey70',stroke =0,size = 1, alpha = 0.6) + 
	geom_smooth(fill=RNcolors[2],color=RNcolors[2],method = 'lm') + 
	geom_text(aes(x=0.04+ min(d$x),y=max(d$y),label = paste('p = ',signif(pdat$nBackP[transInd],2),sep='')),size=2.5) +
	xlab(paste('P(',transLabels[transInd],')',sep='')) +
	ggtitle('n-back') + ylab('WM Performance') + theme_classic() + theme(text = element_text(size = 8)) +
	theme(plot.title = element_text(size = 8,hjust = 0.5))
ggsave(plot = p,filename = paste(savedir,'nBack',transLabels[transInd],'Dprime_k',numClusters,'.pdf',sep =''),units = 'in',height = 2,width = 2.5)
