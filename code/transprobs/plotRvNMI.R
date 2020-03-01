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

clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- getClusterColors(numClusters)
RNcolors <- c('#005C9F','#FF8400')  

# compare transition vs. transition + persistence MI

restMI <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombMutualInfo_k',numClusters,name_root,'.mat',sep = ''))$subjectMI
nbackMI <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombMutualInfo_k',numClusters,name_root,'.mat',sep = ''))$subjectMI
TR <- 3

grps <- rbind(matrix('Rest',nrow = nrow(restMI),ncol = ncol(restMI)),matrix('n-back',nrow = nrow(nbackMI),ncol = ncol(nbackMI)))
p <- ggplot() + geom_violin(aes(x = as.vector(grps), y = as.vector(rbind(restMI,nbackMI)),fill = as.vector(grps))) + theme_classic() +
  scale_fill_manual(limits = c('Rest','n-back'), values = RNcolors) + 
  scale_x_discrete(limits = c('Rest','n-back'),breaks = c('Rest','n-back')) +
  scale_y_continuous(limits = c(0,0.42),breaks = c(0,0.2,0.4)) +
  ylab("Normalized Mutual Information") + xlab("") +theme(text = element_text(size = 8)) +
  theme(legend.position = 'none') + ggtitle('Persistence + Transitions') + theme(plot.title = element_text(size=8,hjust = 0.5))

diffs.full <- t.test(restMI,nbackMI,paired = TRUE)
print(diffs.full)
diffs <- diffs.full$p.value*2 # correct over 2 comparisons for persist + trans and transitions only

print(diffs)
if(diffs < 10^-15){
    p <- p + annotate("text", x = "n-back", y = 1.1*max(rbind(restMI,nbackMI)),label = "**",color = 'red')
}

ggsave(plot = p, filename = paste(masterdir,'analyses/transitionprobabilities/RvNFullMutualInformation3sLag_k',numClusters,'.pdf',sep =""),height = 2,width = 1.8, units = "in")

restMI <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombMutualInfo_k',numClusters,name_root,'.mat',sep = ''))$subjectMITransition
nbackMI <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombMutualInfo_k',numClusters,name_root,'.mat',sep = ''))$subjectMITransition
maxGap <- 1
TR <- 3

grps <- rbind(matrix('Rest',nrow = nrow(restMI),ncol = ncol(restMI)),matrix('n-back',nrow = nrow(nbackMI),ncol = ncol(nbackMI)))
p <- ggplot() + geom_violin(aes(x = as.vector(grps), y = as.vector(rbind(restMI,nbackMI)),fill = as.vector(grps))) + theme_classic() +
  scale_fill_manual(limits = c('Rest','n-back'), values = RNcolors) + 
  scale_x_discrete(limits = c('Rest','n-back'),breaks = c('Rest','n-back')) +
  scale_y_continuous(limits = c(0,0.42),breaks = c(0,0.2,0.4)) +
  ylab("Normalized Mutual Information") + xlab("") +theme(text = element_text(size = 8)) +
  theme(legend.position = 'none') + ggtitle('Transitions Only') + theme(plot.title = element_text(size=8,hjust = 0.5))

diffs.full <- t.test(restMI,nbackMI,paired = TRUE)
print(diffs.full)
diffs <- diffs.full$p.value*2 # correct over 2 comparisons for persist + trans and transitions only
print(diffs)
if(diffs < 10^-15){
    p <- p + annotate("text", x = "n-back", y = 1.1*max(rbind(restMI,nbackMI)),label = "**",color = 'red')
}
ggsave(plot = p, filename = paste(masterdir,'analyses/transitionprobabilities/RvNTransitionMutualInformation3sLag_k',numClusters,'.pdf',sep =""),height = 2,width = 1.8, units = "in")
