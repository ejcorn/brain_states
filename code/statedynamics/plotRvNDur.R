args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))

restDur <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombFractionalOccupancy_k',numClusters,name_root,'.mat',sep = ''))$FractionalOccupancy * 100
nbackDur <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombFractionalOccupancy_k',numClusters,name_root,'.mat',sep = ''))$FractionalOccupancy * 100
clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- getClusterColors(numClusters)
RNcolors <- c('#005C9F','#FF8400')  

grps <- rbind(matrix('Rest',nrow = nrow(restDur),ncol = ncol(restDur)),matrix('n-back',nrow = nrow(nbackDur),ncol = ncol(nbackDur)))
states <- sapply(1:numClusters, function(K) rep(as.character(K),nrow(restDur)))
df.plt <- data.frame(states=as.vector(rbind(states,states)),dur=as.vector(rbind(restDur,nbackDur)),grps=as.vector(grps))
p <- ggplot(df.plt) + geom_split_violin(aes(x = states ,y = dur,fill = grps)) + theme_classic() +
  scale_fill_manual(limits = c('Rest','n-back'), values = RNcolors) + 
  ylab("Fractional Occupancy (%)") + xlab("") +theme(text = element_text(size = 8)) +
  theme(legend.title = element_blank()) + 
  scale_x_discrete(limits = 1:numClusters, breaks=1:numClusters, labels = list(clusterNames)) +
  theme(legend.key.size = unit(0.5,'line'))

diffs.full <- lapply(1:numClusters, function(i) t.test(restDur[,i],nbackDur[,i],paired=TRUE))
print(diffs.full)
diffs <- sapply(diffs.full, function(x) x$p.value)
diffs <- p.adjust(diffs,method = "bonf")
print(diffs)

for(K in 1:numClusters){
  if(diffs[K] < 10^-15){
    p <- p + annotate("text", x = K, y = 1.1*max(rbind(restDur,nbackDur)),label = "**",color = 'red')
  } else if(diffs[K] < 10^-4){
  	p <- p + annotate("text", x = K, y = 1.1*max(rbind(restDur,nbackDur)),label = "*",color = 'red')
  }
}

if(numClusters == 5 | numClusters == 6){
  p <- p + theme(axis.text.x = element_text(size=8,colour = clusterColors))
}
save(df.plt,clusterNames,clusterColors,p,diffs.full,file=paste0(masterdir,'analyses/transitionprobabilities/Fig3a__FractionalOccupancy.RData'))
ggsave(plot = p, filename = paste(masterdir,'analyses/transitionprobabilities/RvNFractionalOccupancy_k',numClusters,'.pdf',sep =""),
  height = 2,width = (numClusters-2) + 0.25, units = "in")