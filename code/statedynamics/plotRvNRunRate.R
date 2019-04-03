args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))

restRunRate <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombDwellTime_k',numClusters,name_root,'.mat',sep = ''))$RunRate
nbackRunRate <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombDwellTime_k',numClusters,name_root,'.mat',sep = ''))$RunRate
clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
RNcolors <- c('#005C9F','#FF8400')  

grps <- rbind(matrix('Rest',nrow = nrow(restRunRate),ncol = ncol(restRunRate)),matrix('n-back',nrow = nrow(nbackRunRate),ncol = ncol(nbackRunRate)))
states <- sapply(1:numClusters, function(K) rep(as.character(K),nrow(restRunRate)))
df <- data.frame(states=as.vector(rbind(states,states)),grps = as.vector(grps),dt=as.vector(rbind(restRunRate,nbackRunRate)))
p <- ggplot(data=df) + geom_split_violin(aes(x = states, y = dt,fill = grps)) + theme_classic() +
  scale_fill_manual(limits = c('Rest','n-back'), values = RNcolors) + 
  ylab("Appearance Rate (runs/min)") + xlab("") +theme(text = element_text(size = 8)) +
  theme(legend.title = element_blank()) + 
  scale_x_discrete(limits = 1:numClusters, breaks=1:numClusters, labels = list(clusterNames)) +
  theme(legend.key.size = unit(0.5,'line'))

diffs.full <- lapply(1:numClusters, function(i) t.test(restRunRate[,i],nbackRunRate[,i],paired=TRUE))
print(diffs.full)
diffs <- sapply(diffs.full, function(x) x$p.value)
diffs <- p.adjust(diffs,method = "bonf")
for(K in 1:numClusters){
  if(diffs[K] < 10^-15){
    p <- p + annotate("text", x = K, y = 1.1*max(rbind(restRunRate,nbackRunRate)),label = "**",color = 'red')
  } else if(diffs[K] < 10^-4){
  	p <- p + annotate("text", x = K, y = 1.1*max(rbind(restRunRate,nbackRunRate)),label = "*",color = 'red')
  }
}

if(numClusters == 5){
  p <- p + theme(axis.text.x = element_text(size=8,colour = clusterColors))
}

ggsave(plot = p, filename = paste(masterdir,'analyses/transitionprobabilities/RvNRunRateTime_k',numClusters,'.pdf',sep =""),height = 2,width = (numClusters-1) + 0.25, units = "in")