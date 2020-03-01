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

savedir <- paste(masterdir,'analyses/development/blockdur/',sep = '')
dir.create(savedir,recursive = TRUE)

################################################
### n-back block specific dwell time and age ###
################################################

# n-back block duration to block-specific dprime

BlockNames <- c('0back','1back','2back')
numBlocks <- length(BlockNames)
restFO <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombFractionalOccupancy_k',
                         numClusters,name_root,'.mat',sep = ''))$FractionalOccupancy
nbackBlockFO <- lapply(1:numBlocks, function(B) readMat(paste(masterdir,'analyses/nbackblocks/FractionalOccupancy',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockFractionalOccupancy)

nbackBlockFO.dprime <- vector("list",numBlocks)	# initialize results list, extract first non-intercept coefficient (reason for index of 2)

for(B in 1:3){
	BlockFO.dprime.output <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBlockFO[[B]][,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
		data = demo)))$coefficients)
	print(paste(B-1,'-back Block',sep=''))
	print(BlockFO.dprime.output)
	nbackBlockFO.dprime[[B]] <- lapply(BlockFO.dprime.output, function(X) X[2,c('Standardized','Pr(>|t|)')])
}
restFO.dprime <- lapply(1:numClusters, function(K) summary(lm.beta(lm(restFO[,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# plot

pdat <- as.data.frame(cbind(t(as.data.frame(nbackBlockFO.dprime[[1]], row.names = c('ZbackB','ZbackP'))),
	t(as.data.frame(nbackBlockFO.dprime[[2]],row.names = c('ObackB','ObackP'))),
	t(as.data.frame(nbackBlockFO.dprime[[3]],row.names = c('TbackB','TbackP'))),
	t(as.data.frame(restFO.dprime,row.names = c('RestB','RestP')))))

p.list <- list.posthoc.correct(list(pdat$ZbackP,pdat$ObackP,pdat$TbackP,pdat$RestP),method = 'bonf')
pdat$ZbackP <- ifelse(p.list[[1]] < 0.05,'*','')
pdat$ObackP <- ifelse(p.list[[2]] < 0.05,'*','')
pdat$TbackP <- ifelse(p.list[[3]] < 0.05,'*','')
pdat$RestP <- ifelse(p.list[[4]] < 0.05,'*','')

rownames(pdat) <- clusterNames
pdat <- data.frame(B = c(pdat$RestB,pdat$ZbackB,pdat$ObackB,pdat$TbackB),
	p=c(pdat$RestP,pdat$ZbackP,pdat$ObackP,pdat$TbackP))
pdat$Scan = rep(c('00Rest','0-back','1-back','2-back'),each=numClusters) # this 00 makes rest appear first in orer of bars
pdat$State = rep(clusterNames,length(p.list))

BlockColors <- c("#FBB4AE", "#B3CDE3", "#CCEBC5") # brewer.pal(3,'Pastel1')
PlotColors <- c(RNcolors[1],BlockColors)
p <- ggplot(data = pdat, aes(y = B, x = State, fill = Scan,label = p)) + 
	geom_bar(stat = "identity", position = position_dodge(width = NULL),alpha = 0.8, color = 'black') +
  geom_text(position = position_dodge(width = 0.9), size = 5) + xlab("") + ylab(expression(beta["Age"])) + ggtitle('Block Fractional Occupancy') +
  scale_fill_manual(values = PlotColors) + scale_x_discrete(limits = clusterNames,breaks = clusterNames) + 
  scale_y_continuous(limits = c(1.05*min(pdat$B),1.1*max(pdat$B)))+ 
  theme_classic() + theme(text = element_text(size = 8)) + theme(legend.key.size = unit(0.5,'line')) + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.95)) +
  theme(legend.title = element_blank(),plot.title = element_text(hjust =0.5,size = 8),legend.position = 'none') +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
if(numClusters == 5){
	p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

save(p,pdat,file=paste0(savedir,'Fig6a__RestBlockFOAge.RData'))
ggsave(plot = p,filename = paste(savedir,'RestAndnBackBlockAgeFO_k',numClusters,'.pdf',sep =''),units = 'cm',height = 3,width = 6)