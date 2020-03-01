args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

library(ggplot2)
library(R.matlab)
library(RColorBrewer)
library(plotrix)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))

BlockNames <- c('0back','1back','2back')
numBlocks <- length(BlockNames)

clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
#clusterNames <- c('A',"b","c","d","e")
clusterColors <- getClusterColors(numClusters)
RNcolors <- c('#005C9F','#FF8400') 

restDwell <- colMeans(readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombDwellTime_k',
	numClusters,name_root,'.mat',sep = ''))$DwellTimeMean)
restSEM <- apply(X = readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombDwellTime_k',
	numClusters,name_root,'.mat',sep = ''))$DwellTimeMean,FUN = std.error, MARGIN= 2)
nbackDwell <- lapply(1:numBlocks, function(B) colMeans(readMat(paste(masterdir,'analyses/nbackblocks/DwellTime',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockDwellTime))
nbackDwellSEM <- lapply(1:numBlocks, function(B) apply(X = readMat(paste(masterdir,'analyses/nbackblocks/DwellTime',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockDwellTime,FUN=std.error,MARGIN=2))
nbackDwell <- c(list(restDwell),nbackDwell)
nbackDwellSEM <- c(list(restSEM),nbackDwellSEM)
numBlocks <- numBlocks + 1

states <- lapply(1:numBlocks, function(B) sapply(1:numClusters, function(K) as.character(K)))
blocks <- lapply(1:numBlocks, function(B) sapply(1:numClusters, function(K) B))
nbackDwell <- as.vector(do.call(cbind,nbackDwell))
nbackDwellSEM <- as.vector(do.call(cbind,nbackDwellSEM))
states <- as.vector(do.call(cbind,states))
blocks <- as.vector(do.call(cbind,blocks))

BlockNames <- c('Rest','0-back','1-back','2-back')

ylim <- c(min(nbackDwell - 2*nbackDwellSEM),max(nbackDwell + 2*nbackDwellSEM))
ylim.stretch <- 0.05*(ylim[2] - ylim[1])*c(-1,1)     # amount by which to pad y limits
ylim <- ylim + ylim.stretch
p <- ggplot() + geom_line(aes(x = blocks, y = nbackDwell, color = states), size = 0.4) +
  geom_errorbar(aes(ymin = nbackDwell - 2*nbackDwellSEM,ymax = nbackDwell + 2*nbackDwellSEM, x = blocks, color = states),width = 0.1,size = 0.25) +
  scale_x_continuous(limits = c(0.7,numBlocks+0.25),breaks= c(1:numBlocks),label=BlockNames,expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),limits=ylim) + 
  xlab("Task Block") + ylab("Dwell Time (seconds)") + 
  theme_classic() + theme(text = element_text(size = 8),legend.position = 'none')

if(numClusters == 5 | numClusters == 6){
  p <- p + scale_color_manual(limits = c(1:numClusters), values = clusterColors, label=clusterNames) +
    annotate("text",x = rep(1,numClusters),y = 0.1 + (nbackDwell + nbackDwellSEM)[blocks == 1],label = clusterNames,size = 2, color = clusterColors)
  
} else {
  #p <- p + annotate("text",x = rep(1,numClusters),y = 0.7 + (nbackDwell + nbackDwellSEM)[blocks == 1],label = clusterNames,size = 2,color = unique(states)) +
  p <- p + geom_text(aes(x = rep(1,numClusters),y = 0.1 + (nbackDwell + nbackDwellSEM)[blocks == 1],label = clusterNames,color = unique(states)),size = 2) +
    scale_color_discrete(limits = c(1:numClusters))
}

ggsave(plot = p, filename = paste(masterdir,'analyses/nbackblocks/nbackBlockDwellTime_k',numClusters,'.pdf',sep =""),
  height = 2,width = 2.25, units = "in")