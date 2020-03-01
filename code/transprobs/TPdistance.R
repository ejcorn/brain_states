args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

library(R.matlab)
library(ggplot2)
library(lm.beta)
masterdir <- paste(basedir,'results/',name_root,'/',sep='')
source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))
scanlab <- c("RestComb","nBackComb")
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep =""))

clusterAssignments <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterAssignments$clusterAssignments[[1]][[5]])
kClusterCentroids <- clusterAssignments$clusterAssignments[[1]][[2]]

savedir = paste(masterdir,'/analyses/control_energy',sep='');
dir.create(savedir,recursive = TRUE)

# account for distance between centroids

persistExclude = (1:numClusters) + (numClusters*(0:(numClusters-1)))
InterStateDistance = matrix(0,nrow = numClusters,ncol = numClusters)
for(K1 in 1:numClusters){
	for(K2 in 1:numClusters){
		InterStateDistance[K1,K2] = sqrt(sum((kClusterCentroids[,K1] - kClusterCentroids[,K2])^2))
	}
}

InterStateDistance = as.vector(t(InterStateDistance))[-persistExclude]
rtp <- colMeans(readMat(paste(basedir,"results/",name_root,"/analyses/transitionprobabilities/",
	scanlab[1],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability)[-persistExclude]
cor.test(InterStateDistance,rtp)

ntp <- colMeans(readMat(paste(basedir,"results/",name_root,"/analyses/transitionprobabilities/",
	scanlab[2],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability)[-persistExclude]
cor.test(InterStateDistance,ntp)

rtp <- readMat(paste(basedir,"results/",name_root,"/analyses/transitionprobabilities/",
	scanlab[1],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability[,-persistExclude]

ntp <- readMat(paste(basedir,"results/",name_root,"/analyses/transitionprobabilities/",
	scanlab[2],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability[,-persistExclude]

nobs <- nrow(rtp)
rest.distcors <- lapply(1:nobs, function(N) cor.test(InterStateDistance,rtp[N,]))
rest.rp <- sapply(rest.distcors,function(d) (d$p.value < 0.05)*d$estimate)
rest.r <- sapply(rest.distcors,function(d) d$estimate)
rest.rsig <- sapply(rest.distcors,function(d) (d$p.value < 0.05))
nback.distcors <- lapply(1:nobs, function(N) cor.test(InterStateDistance,ntp[N,]))
nback.rp <- sapply(nback.distcors,function(d) (d$p.value < 0.05)*d$estimate)
nback.r <- sapply(nback.distcors,function(d) d$estimate)
nback.rsig <- sapply(nback.distcors,function(d) (d$p.value < 0.05))

dat <- c(rest.r,nback.r)
grps <- c(rep('Rest',nobs),rep('n-back',nobs))
df <- data.frame(dat=dat,grps=grps)

p <- ggplot(df) + geom_violin(aes(x=grps,y=dat,fill=grps)) + theme_classic() +
  scale_fill_manual(limits = c('Rest','n-back'), values = RNcolors) + 
  ylab("r(TP,Distance)") + xlab("") +theme(text = element_text(size = 8)) + 
  scale_x_discrete(limits = c('Rest','n-back'),breaks= c('Rest','n-back'),label =c('Rest','n-back')) +
  theme(legend.position = 'none') + theme(axis.ticks.x = element_blank())
p

test.output <- t.test(nback.r,rest.r,paired = TRUE)
print(test.output)
if(test.output$p.value < 10^-15){
	p <- p + annotate("text", x = 'n-back', y = 0.5,label = "**",color = 'red')
}

ggsave(plot = p,filename = paste(savedir,'/RvNDistanceTPCorr_k',numClusters,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2)