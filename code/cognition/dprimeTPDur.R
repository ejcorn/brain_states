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

nbbeh <- read.csv('/data/jag/bassett-lab/Lausanne1601/nback/n1601_nbackBehavior_from_20160207_dataRelease.csv')
nbbeh <- nbbeh[which(nbbeh$scanid %in% demo$scanid),]
nbbeh <- nbbeh[order(nbbeh$scanid),]

cnb <- read.csv('/data/jag/bassett-lab/Lausanne1601/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv')
cnb <- cnb[which(nbbeh$scanid %in% demo$scanid),]
cnb <- cnb[order(cnb$scanid),]

data <- cbind(demo, cnb[,grepl('Accuracy',colnames(cnb))], nbbeh[,grepl('Dprime',colnames(nbbeh))])

#################################
### Cognition and Dwell Times ###
#################################

savedir <- paste(masterdir,'analyses/cognition/dwelltime/',sep = '')
if(!dir.exists(savedir)){
	dir.create(savedir,recursive = TRUE)
}

restDur <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombStateDurations_k',
                         numClusters,name_root,'.mat',sep = ''))$stateDuration * 100
nbackREDur <- readMat(paste(masterdir,'analyses/nbackblocks/nBackRestExcludeDwellTime_k',
                          numClusters,name_root,'.mat',sep = ''))$BlockDwellTime*100
#nbackREDur <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombStateDurations_k',
#                        numClusters,name_root,'.mat',sep = ''))$stateDuration * 100

# rest and n-back to overall d-prime

restDur.dprime <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBehAllDprime ~ restDur[,K] + age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackREDur.dprime <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBehAllDprime ~ nbackREDur[,K] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# plot

pdat <- as.data.frame(cbind(t(as.data.frame(restDur.dprime, row.names = c('RestB','RestP'))),t(as.data.frame(nbackREDur.dprime,row.names = c('nBackB','nBackP')))))
pdat$RestP <- p.adjust(pdat$RestP,method = 'bonf') < 0.05
pdat$nBackP <- p.adjust(pdat$nBackP,method = 'bonf') < 0.05
pdat$RestP <- ifelse(pdat$RestP,'*','')
pdat$nBackP <- ifelse(pdat$nBackP,'*','')
colnames(pdat) <- rep(c('B','p'),2)
pdat <- as.data.frame(rbind(pdat[1:numClusters,1:2],pdat[1:numClusters,3:4]))
pdat$Scan = c(rep('Rest',numClusters),rep('n-back',numClusters))
pdat$State = rep(clusterNames,2)

p <- ggplot(data = pdat, aes(y = B, x = State, fill = Scan,label = p)) + 
	geom_bar(stat = "identity", position = position_dodge(width = NULL), color = 'black') +
  geom_text(position = position_dodge(width = 0.9), size = 5) + xlab("") + ylab(expression(beta["DT"])) + ggtitle('Overall WM Performance') +
  scale_y_continuous(limits = c(1.05*min(pdat$B),1.1*max(pdat$B))) +
  scale_fill_manual(limits = c('Rest','n-back'), values = RNcolors) + scale_x_discrete(limits = clusterNames,breaks = clusterNames) + 
  theme_classic() + theme(text = element_text(size = 8)) + theme(legend.position = 'none',plot.title = element_text(hjust =0.5,size = 8)) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.95)) + theme(legend.key.size = unit(0.5,'line')) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

if(numClusters == 5){
	p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

ggsave(plot = p,filename = paste(savedir,'RvNOverallDprimeDwellTime_k',numClusters,'.pdf',sep =''),units = 'cm',height = 3,width = 5)

# rest and n-back to overall accuracy

restDur.acc <- lapply(1:numClusters, function(K) summary(lm.beta(lm(Overall_Accuracy ~ restDur[,K] + age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackREDur.acc <- lapply(1:numClusters, function(K) summary(lm.beta(lm(Overall_Accuracy ~ nbackREDur[,K] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# plot

pdat <- as.data.frame(cbind(t(as.data.frame(restDur.acc, row.names = c('RestB','RestP'))),t(as.data.frame(nbackREDur.acc,row.names = c('nBackB','nBackP')))))
pdat$RestP <- p.adjust(pdat$RestP,method = 'bonf') < 0.05
pdat$nBackP <- p.adjust(pdat$nBackP,method = 'bonf') < 0.05
pdat$RestP <- ifelse(pdat$RestP,'*','')
pdat$nBackP <- ifelse(pdat$nBackP,'*','')
colnames(pdat) <- rep(c('B','p'),2)
pdat <- as.data.frame(rbind(pdat[1:numClusters,1:2],pdat[1:numClusters,3:4]))
pdat$Scan = c(rep('Rest',numClusters),rep('n-back',numClusters))
pdat$State = rep(clusterNames,2)

p <- ggplot(data = pdat, aes(y = B, x = State, fill = Scan,label = p)) + 
	geom_bar(stat = "identity", position = position_dodge(width = NULL), color = 'black') +
  geom_text(position = position_dodge(width = 0.9), size = 5) + xlab("") + ylab(expression(beta["DT"])) + ggtitle('Overall Accuracy') +
  scale_fill_manual(limits = c('Rest','n-back'), values = RNcolors) + scale_x_discrete(limits = clusterNames,breaks = clusterNames) + 
  scale_y_continuous(limits = c(1.05*min(pdat$B),1.1*max(pdat$B))) + theme(legend.position = 'none') +
  theme_classic() + theme(text = element_text(size = 8)) + theme(legend.title = element_blank(),plot.title = element_text(hjust =0.5)) + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

if(numClusters == 5){
	p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

ggsave(plot = p,filename = paste(savedir,'RvNOverallAccuracyDwellTime_k',numClusters,'.pdf',sep =''),units = 'cm',height = 3,width = 5)

# n-back block duration to block-specific dprime

BlockNames <- c('0back','1back','2back')
numBlocks <- length(BlockNames)
nbackBlockDur <- lapply(1:numBlocks, function(B) readMat(paste(masterdir,'analyses/nbackblocks/DwellTime',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockDwellTime*100)

nbackBlockDur.dprime <- vector("list",numBlocks)
nbackBlockDur.dprime[[1]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBehZerobackDprime ~ nbackBlockDur[[1]][,K] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackBlockDur.dprime[[2]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBehOnebackDprime ~ nbackBlockDur[[2]][,K] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackBlockDur.dprime[[3]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBehTwobackDprime ~ nbackBlockDur[[3]][,K] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# plot

pdat <- as.data.frame(cbind(t(as.data.frame(nbackBlockDur.dprime[[1]], row.names = c('ZbackB','ZbackP'))),
	t(as.data.frame(nbackBlockDur.dprime[[2]],row.names = c('ObackB','ObackP'))),
	t(as.data.frame(nbackBlockDur.dprime[[3]],row.names = c('TbackB','TbackP')))))
pdat$ZbackP <- p.adjust(pdat$ZbackP,method = 'bonf') < 0.05
pdat$ObackP <- p.adjust(pdat$ObackP,method = 'bonf') < 0.05
pdat$TbackP <- p.adjust(pdat$TbackP,method = 'bonf') < 0.05
pdat$ZbackP <- ifelse(pdat$ZbackP,'*','')
pdat$ObackP <- ifelse(pdat$ObackP,'*','')
pdat$TbackP <- ifelse(pdat$TbackP,'*','')
colnames(pdat) <- rep(c('B','p'),3)
pdat <- as.data.frame(rbind(pdat[1:numClusters,1:2],pdat[1:numClusters,3:4],pdat[1:numClusters,5:6]))
pdat$Scan = c(rep('0-back',numClusters),rep('1-back',numClusters),rep('2-back',numClusters))
pdat$State = rep(clusterNames,3)

p <- ggplot(data = pdat, aes(y = B, x = State, fill = Scan,label = p)) + 
	geom_bar(stat = "identity", position = position_dodge(width = NULL),alpha = 0.8, color = 'black') +
  geom_text(position = position_dodge(width = 0.9), size = 5) + xlab("") + ylab(expression(beta["DT"])) + ggtitle('Block WM Performance') +
  scale_fill_brewer(limits = unique(pdat$Scan), palette = 'Pastel1') + scale_x_discrete(limits = clusterNames,breaks = clusterNames) + 
  scale_y_continuous(limits = c(1.05*min(pdat$B),1.1*max(pdat$B)))+ 
  theme_classic() + theme(text = element_text(size = 8)) + theme(legend.key.size = unit(0.5,'line')) + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.95)) +
  theme(legend.title = element_blank(),plot.title = element_text(hjust =0.5,size = 8),legend.position = 'none') +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
if(numClusters == 5){
	p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

ggsave(plot = p,filename = paste(savedir,'RvNBlockDprimeDwellTime_k',numClusters,'.pdf',sep =''),units = 'cm',height = 3,width = 6)

##############################
### Cognition and Dynamics ###
##############################

symmvars <- readMat(paste(masterdir,'/analyses/transitionprobabilities/symmetry/RvNSymmetryScorev2_k',
		numClusters,'.mat',sep = ''))
restSymm <- symmvars$restSymmetryScore
nbackSymm <-symmvars$nBackSymmetryScore

summary(lm(nbackBehAllDprime ~ restSymm + age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex,
	data = data))
summary(lm(nbackBehAllDprime ~ nbackSymm + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = data))

# MI <- readMat(paste(masterdir,'analyses/transitionprobabilities/MutualInfoGap1to6_k',numClusters,name_root,'.mat',sep = ''))$subjectMI
# restMI <- as.matrix(MI[[1]][[1]][,1])
# nbackMI <- as.matrix(MI[[2]][[1]][,1])

# summary(lm(nbackBehAllDprime ~ restMI + age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex,
# 	data = data))
# summary(lm(nbackBehAllDprime ~ nbackMI + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
# 	data = data))

#################################
### Cognition and Trans Probs ###
#################################

savedir <- paste(masterdir,'analyses/cognition/transprobs/',sep = '')
if(!dir.exists(savedir)){
	dir.create(path = savedir,recursive = TRUE)
}

restTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombTransitionProbabilities_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability
nbackTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombTransitionProbabilities_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability

numTransitions <- ncol(restTP)
transLabels <- as.vector(sapply(1:numClusters, function(i) sapply(1:numClusters,function(j) paste(clusterNames[i],'to',clusterNames[j]))))	#label transitions

# overall d prime
restTP.dprime <- lapply(1:numTransitions, function(T) summary(lm.beta(lm(nbackBehAllDprime ~ restTP[,T] + age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackTP.dprime <- lapply(1:numTransitions, function(T) summary(lm.beta(lm(nbackBehAllDprime ~ nbackTP[,T] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])

pdat <- as.data.frame(t(as.data.frame(restTP.dprime, row.names = c('RestB','RestP'))))
pdat$RestP <- p.adjust(pdat$RestP,method = 'bonf') < 0.05
v <- matrix((pdat$RestB*pdat$RestP), nrow = numClusters, byrow = TRUE)
p <- TP.beta.plot(v,clusterColors,title = 'Rest: d-prime')
ggsave(plot = p,filename = paste(savedir,'RestTPDprime_k',numClusters,'.pdf',sep =''),units = 'cm',height = 4,width = 4)

pdat <- as.data.frame(t(as.data.frame(nbackTP.dprime, row.names = c('nBackB','nBackP'))))
pdat$nBackP <- p.adjust(pdat$nBackP,method = 'fdr') < 0.05
v <- matrix((pdat$nBackB*pdat$nBackP), nrow = numClusters, byrow = TRUE)
p <- TP.beta.plot(v,clusterColors,title = 'n-back: d-prime')
ggsave(plot = p,filename = paste(savedir,'nBackTPDprime_k',numClusters,'.pdf',sep =''),units = 'cm',height = 4,width = 4)


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


