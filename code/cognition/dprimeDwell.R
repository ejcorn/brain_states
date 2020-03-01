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

################################
### Cognition and Dwell Time ###
################################

savedir <- paste(masterdir,'analyses/cognition/dwelltime/',sep = '')
dir.create(savedir,recursive = TRUE)

restDwell <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombDwellTime_k',
                         numClusters,name_root,'.mat',sep = ''))$DwellTimeMean

# rest to overall d-prime

restDwell.dprime <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBehAllDprime ~ restDwell[,K] + age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# plot

pdat <- as.data.frame(t(as.data.frame(restDwell.dprime, row.names = c('RestB','RestP'))))
pdat$RestP <- p.adjust(pdat$RestP,method = 'bonf') < 0.05
pdat$RestP <- ifelse(pdat$RestP,'*','')
colnames(pdat) <- rep(c('B','p'),1)
pdat$Scan = rep('Rest',numClusters)
pdat$State = clusterNames

p <- ggplot(data = pdat, aes(y = B, x = State, fill = Scan,label = p)) + 
	geom_bar(stat = "identity", position = position_dodge(width = NULL), color = 'black') +
  geom_text(position = position_dodge(width = 0.9), size = 5) + xlab("") + ylab(expression(beta["DT"])) + ggtitle('Overall WM Performance') +
  scale_y_continuous(limits = c(1.05*min(pdat$B),1.1*max(pdat$B))) +
  scale_fill_manual(limits = c('Rest'), values = RNcolors[1]) + scale_x_discrete(limits = clusterNames,breaks = clusterNames) + 
  theme_classic() + theme(text = element_text(size = 8)) + theme(legend.position = 'none',plot.title = element_text(hjust =0.5,size = 8)) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.95)) + theme(legend.key.size = unit(0.5,'line')) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

if(numClusters == 5){
	p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

ggsave(plot = p,filename = paste(savedir,'RestOverallDprimeDwellTime_k',numClusters,'.pdf',sep =''),units = 'cm',height = 3,width = 5)

# n-back block duration to block-specific dprime

BlockNames <- c('0back','1back','2back')
numBlocks <- length(BlockNames)
nbackBlockDwell <- lapply(1:numBlocks, function(B) readMat(paste(masterdir,'analyses/nbackblocks/DwellTime',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockDwellTime)

nbackBlockDwell.dprime <- vector("list",numBlocks)
nbackBlockDwell.dprime[[1]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBehZerobackDprime ~ nbackBlockDwell[[1]][,K] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackBlockDwell.dprime[[2]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBehOnebackDprime ~ nbackBlockDwell[[2]][,K] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackBlockDwell.dprime[[3]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBehTwobackDprime ~ nbackBlockDwell[[3]][,K] + age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = data)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# plot

pdat <- as.data.frame(cbind(t(as.data.frame(nbackBlockDwell.dprime[[1]], row.names = c('ZbackB','ZbackP'))),
	t(as.data.frame(nbackBlockDwell.dprime[[2]],row.names = c('ObackB','ObackP'))),
	t(as.data.frame(nbackBlockDwell.dprime[[3]],row.names = c('TbackB','TbackP')))))
# pdat$ZbackP <- p.adjust(pdat$ZbackP,method = 'bonf') < 0.05
# pdat$ObackP <- p.adjust(pdat$ObackP,method = 'bonf') < 0.05
# pdat$TbackP <- p.adjust(pdat$TbackP,method = 'bonf') < 0.05
# pdat$ZbackP <- ifelse(pdat$ZbackP,'*','')
# pdat$ObackP <- ifelse(pdat$ObackP,'*','')
# pdat$TbackP <- ifelse(pdat$TbackP,'*','')
p.list <- list.posthoc.correct(list(pdat$ZbackP,pdat$ObackP,pdat$TbackP),method = 'bonf')
pdat$ZbackP <- ifelse(p.list[[1]] < 0.05,'*','')
pdat$ObackP <- ifelse(p.list[[2]] < 0.05,'*','')
pdat$TbackP <- ifelse(p.list[[3]] < 0.05,'*','')

colnames(pdat) <- rep(c('B','p'),3)
pdat <- as.data.frame(rbind(pdat[1:numClusters,1:2],pdat[1:numClusters,3:4],pdat[1:numClusters,5:6]))
pdat$Scan = c(rep('0-back',numClusters),rep('1-back',numClusters),rep('2-back',numClusters))
pdat$State = rep(clusterNames,3)

p <- ggplot(data = pdat, aes(y = B, x = State, fill = Scan,label = p)) + 
	geom_bar(stat = "identity", position = position_dodge(width = NULL),alpha = 0.8, color = 'black',size=0.5) +
  geom_text(position = position_dodge(width = 0.9), size = 5) + xlab("") + ylab(expression(beta["DT"])) + ggtitle('Block WM Performance') +
  scale_fill_brewer(limits = unique(pdat$Scan), palette = 'Pastel1') + scale_x_discrete(limits = clusterNames,breaks = clusterNames) + 
  scale_y_continuous(limits = c(1.05*min(pdat$B),1.1*max(pdat$B)))+ 
  theme_classic() + theme(text = element_text(size = 8)) + theme(legend.key.size = unit(0.5,'line')) + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.95)) +
  theme(legend.title = element_blank(),plot.title = element_text(hjust =0.5,size = 8),legend.position = 'none') +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

height = 3; width = 6
if(numClusters == 5 | numClusters == 6){
  	p <- p + theme(axis.text.y = element_text(color = rev(clusterColors))) 
    p <- p + coord_flip()
    p <- p + scale_x_discrete(limits = rev(clusterNames),breaks = rev(clusterNames))
    height = 5; width = 2.5
}

ggsave(plot = p,filename = paste(savedir,'RvNBlockDprimeDwellTime_k',numClusters,'.pdf',sep =''),units = 'cm',height = height,width = width)
