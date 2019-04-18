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

savedir <- paste(masterdir,'analyses/development/blockdur/',sep = '')
dir.create(savedir,recursive = TRUE)

################################################
### n-back block specific dwell time and age ###
################################################

# n-back block duration to block-specific dprime

BlockNames <- c('0back','1back','2back')
numBlocks <- length(BlockNames)
nbackBlockDur <- lapply(1:numBlocks, function(B) readMat(paste(masterdir,'analyses/nbackblocks/FractionalOccupancy',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockFractionalOccupancy)

nbackBlockDur.dprime <- vector("list",numBlocks)	# initialize results list, extract first non-intercept coefficient (reason for index of 2)
nbackBlockDur.dprime[[1]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBlockDur[[1]][,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackBlockDur.dprime[[2]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBlockDur[[2]][,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackBlockDur.dprime[[3]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBlockDur[[3]][,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# plot

pdat <- as.data.frame(cbind(t(as.data.frame(nbackBlockDur.dprime[[1]], row.names = c('ZbackB','ZbackP'))),
	t(as.data.frame(nbackBlockDur.dprime[[2]],row.names = c('ObackB','ObackP'))),
	t(as.data.frame(nbackBlockDur.dprime[[3]],row.names = c('TbackB','TbackP')))))
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

if(name_root == 'ScanCLaus250Z0final' && numClusters == 5){
	save(pdat, file=paste(savedir,'Fig7d__nBackBlockAgeFractionalOccupancy_k',numClusters,'.RData',sep =''))
}

p <- ggplot(data = pdat, aes(y = B, x = State, fill = Scan,label = p)) + 
	geom_bar(stat = "identity", position = position_dodge(width = NULL),alpha = 0.8, color = 'black') +
  geom_text(position = position_dodge(width = 0.9), size = 5) + xlab("") + ylab(expression(beta["Age"])) + ggtitle('Block Fractional Occupancy') +
  scale_fill_brewer(limits = unique(pdat$Scan), palette = 'Pastel1') + scale_x_discrete(limits = clusterNames,breaks = clusterNames) + 
  scale_y_continuous(limits = c(1.05*min(pdat$B),1.1*max(pdat$B)))+ 
  theme_classic() + theme(text = element_text(size = 8)) + theme(legend.key.size = unit(0.5,'line')) + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.95)) +
  theme(legend.title = element_blank(),plot.title = element_text(hjust =0.5,size = 8),legend.position = 'none') +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
if(numClusters == 5){
	p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

ggsave(plot = p,filename = paste(savedir,'nBackBlockAgeFractionalOccupancy_k',numClusters,'.pdf',sep =''),units = 'cm',height = 3,width = 6)