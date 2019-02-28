args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])

library(ggplot2)
library(R.matlab)
library(RColorBrewer)
library(plotrix)

#name_root <- 'ScanCLaus250Z2sqeuc'
#numClusters <- 5

masterdir <- paste("/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_",name_root,"/",sep="")
#masterdir <- paste("~/Dropbox/Cornblath_Bassett_Projects/code/control_fc/restnbackpipeline/results/",name_root,'/',sep='')

BlockNames <- c('0back','1back','2back','Rest')
numBlocks <- length(BlockNames)

clusterNames <- readMat(paste('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
#clusterNames <- c('A',"b","c","d","e")
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
RNcolors <- c('#005C9F','#FF8400') 

restDur <- colMeans(readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombStateDurations_k',
	numClusters,name_root,'.mat',sep = ''))$stateDuration * 100)
restSEM <- apply(X = readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombStateDurations_k',
	numClusters,name_root,'.mat',sep = ''))$stateDuration * 100,FUN = std.error, MARGIN= 2)
nbackDur <- lapply(1:numBlocks, function(B) colMeans(readMat(paste(masterdir,'analyses/nbackblocks/DwellTime',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockDwellTime*100))
nbackDurSEM <- lapply(1:numBlocks, function(B) apply(X = readMat(paste(masterdir,'analyses/nbackblocks/DwellTime',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockDwellTime*100,FUN=std.error,MARGIN=2))
nbackDur <- c(list(restDur),nbackDur)
nbackDurSEM <- c(list(restSEM),nbackDurSEM)
numBlocks <- numBlocks + 1

states <- lapply(1:numBlocks, function(B) sapply(1:numClusters, function(K) as.character(K)))
blocks <- lapply(1:numBlocks, function(B) sapply(1:numClusters, function(K) B))
nbackDur <- as.vector(do.call(cbind,nbackDur))
nbackDurSEM <- as.vector(do.call(cbind,nbackDurSEM))
states <- as.vector(do.call(cbind,states))
blocks <- as.vector(do.call(cbind,blocks))

BlockNames <- c('Rest','0-back','1-back','2-back','Break')

p <- ggplot() + geom_line(aes(x = blocks, y = nbackDur, color = states), size = 0.4) +
  geom_errorbar(aes(ymin = nbackDur - 2*nbackDurSEM,ymax = nbackDur + 2*nbackDurSEM, x = blocks, color = states),width = 0.1,size = 0.25) +
  scale_x_continuous(limits = c(0.7,numBlocks+0.25),breaks= c(1:numBlocks),label=BlockNames,expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) + 
  xlab("Task Block") + ylab("Time Spent (%)") + 
  theme_classic() + theme(text = element_text(size = 8),legend.position = 'none')

if(numClusters == 5){
  p <- p + scale_color_manual(limits = c(1:numClusters), values = clusterColors, label=clusterNames) +
    annotate("text",x = rep(1,numClusters),y = 0.7 + (nbackDur + nbackDurSEM)[blocks == 1],label = clusterNames,size = 2, color = clusterColors)
  
} else {
  #p <- p + annotate("text",x = rep(1,numClusters),y = 0.7 + (nbackDur + nbackDurSEM)[blocks == 1],label = clusterNames,size = 2,color = unique(states)) +
  p <- p + geom_text(aes(x = rep(1,numClusters),y = 0.7 + (nbackDur + nbackDurSEM)[blocks == 1],label = clusterNames,color = unique(states)),size = 2) +
    scale_color_discrete(limits = c(1:numClusters))
}

ggsave(plot = p, filename = paste(masterdir,'analyses/nbackblocks/nbackBlockTimeSpent_k',numClusters,'.pdf',sep =""),height = 2,width = 2.7, units = "in")