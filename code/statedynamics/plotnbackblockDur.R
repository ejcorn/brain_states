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

colMedians <- function(x){
  return(apply(x,2, function(y) median(y,na.rm = T)))
}

# prepare data for boxplot
restDur <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombFractionalOccupancy_k',
  numClusters,name_root,'.mat',sep = ''))$FractionalOccupancy * 100
restDur <- cbind(matrix.to.df(restDur,dnames=list(NULL,as.character(1:numClusters))),grp=1)
nbackDur <- lapply(1:numBlocks, function(B) readMat(paste(masterdir,'analyses/nbackblocks/FractionalOccupancy',
  BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockFractionalOccupancy*100)
nbackDur <- lapply(1:numBlocks, function(B) cbind(matrix.to.df(nbackDur[[B]],dnames=list(NULL,as.character(1:numClusters))),grp=B+1))
df.FO <- rbind(restDur,do.call('rbind',nbackDur))
df.FO <- collapse.columns(df.FO,cnames=as.character(1:numClusters),groupby='grp')
df.FO$names <- as.character(df.FO$names)

# prepare data for line plot
restDur.mean <- colMeans(readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombFractionalOccupancy_k',
	numClusters,name_root,'.mat',sep = ''))$FractionalOccupancy * 100)
restSEM <- apply(X = readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombFractionalOccupancy_k',
	numClusters,name_root,'.mat',sep = ''))$FractionalOccupancy * 100,FUN = std.error, MARGIN= 2)
nbackDur.mean <- lapply(1:numBlocks, function(B) colMeans(readMat(paste(masterdir,'analyses/nbackblocks/FractionalOccupancy',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockFractionalOccupancy*100))
nbackDurSEM <- lapply(1:numBlocks, function(B) apply(X = readMat(paste(masterdir,'analyses/nbackblocks/FractionalOccupancy',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockFractionalOccupancy*100,FUN=std.error,MARGIN=2))
nbackDur.mean <- c(list(restDur.mean),nbackDur.mean)
nbackDurSEM <- c(list(restSEM),nbackDurSEM)
numBlocks <- numBlocks + 1

states <- lapply(1:numBlocks, function(B) sapply(1:numClusters, function(K) as.character(K)))
blocks <- lapply(1:numBlocks, function(B) sapply(1:numClusters, function(K) B))
nbackDur.mean <- as.vector(do.call(cbind,nbackDur.mean))
nbackDurSEM <- as.vector(do.call(cbind,nbackDurSEM))
states <- as.vector(do.call(cbind,states))
blocks <- as.vector(do.call(cbind,blocks))

BlockNames <- c('Rest','0-back','1-back','2-back')

ylim <- c(min(nbackDur.mean - 2*nbackDurSEM),max(nbackDur.mean + 2*nbackDurSEM))
ylim.stretch <- 0.05*(ylim[2] - ylim[1])*c(-1,1)     # amount by which to pad y limits
ylim <- ylim + ylim.stretch
df.plt <- data.frame(blks=blocks,dur=nbackDur.mean,names=states)
df.FO$group <- as.factor(df.FO$group)
df.FO$names <- as.factor(df.FO$names)
cluster.labs <- clusterNames
names(cluster.labs) <- levels(df.FO$names)
p <- ggplot() + geom_boxplot(data=df.FO,aes(x = group, y = values, fill=names),outlier.shape=NA,fatten=0.7,lwd=0.3,width=0.5) + 
  geom_line(data=df.plt,aes(x=blks,y=dur,color=names),color='black',linetype='dashed',size=0.2)+
  facet_grid(names~.,labeller = labeller(names = cluster.labs),scales='free') +
  scale_x_discrete(limits=levels(df.FO$group),breaks= levels(df.FO$group),label=BlockNames) +
  scale_y_continuous(expand=c(0,0)) +   
  scale_fill_manual(limits = levels(df.FO$names), values = clusterColors, label=clusterNames) +
  xlab("Task Block") + ylab("Fractional Occupancy (%)") + 
  theme_classic() + theme(text = element_text(size = 8),legend.position = 'none',
    strip.text=element_text(size=6,color=clusterColors))
save(df.plt,df.FO,p,file=paste0(masterdir,'analyses/nbackblocks/Fig3d__nbackblocks.RData'))
ggsave(plot = p, filename = paste(masterdir,'analyses/nbackblocks/nbackBlockFractionalOccupancy_k',numClusters,'.pdf',sep =""),
  height = 2,width = 2, units = "in")