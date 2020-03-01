args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]
c <- as.numeric(args[4])
T <- as.numeric(args[5])

# variables that get passed in:
# c: normalization factor
# normalize each matrix by 1/c+max(eig(A)), then subtract identity. ensures max eigenvalue is [0,-1), such that matrices are stable and 
# at long time horizons all nodes either go to 0 or to a constant non-zero eigenvector associated with max eig
# T: control horizon, time over which to exert control

masterdir <- paste(basedir,'results/',name_root,'/',sep='')

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))
source(paste(basedir,'code/miscfxns/statfxns.R',sep=''))

scanlab <- c("RestComb","nBackComb")
persistExclude = (1:numClusters) + (numClusters*(0:(numClusters-1)))
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  
clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- getClusterColors(numClusters)

savedir = paste(masterdir,'/analyses/control_energy',sep='');

# compare energy of maintaining each state with sphere permuted states for actual network and all null models
# .Null indicates sphere permuted null states, .mio indicates null matrix obtained through randmio_und.m, .DLW = Rick Betzel geometric + topological null
matResults <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_FA_k',numClusters,'c',c,'T',T,'.mat',sep = ''))
energy <- matResults$Epersist
null_energy <- matResults$Epersist.Null
energy.randmio <- matResults$Epersist.mio   
null_energy.randmio <- matResults$Epersist.Null.mio
energy.DLW <- matResults$Epersist.DLW
null_energy.DLW  <- matResults$Epersist.Null.DLW

# make vectors to plot
nperms <- nrow(null_energy)
state <- rep(as.vector(sapply(1:numClusters, function(K) rep(as.character(K),nperms))),3)
E <- c(as.vector(null_energy),as.vector(null_energy.randmio),as.vector(null_energy.DLW))
grp <- factor(c(rep('SC',nperms*numClusters),rep('Deg. Pres.',nperms*numClusters),rep('Space Pres.',nperms*numClusters)),levels = c('Deg. Pres.','Space Pres.','SC'),ordered = TRUE)
E.actual <- c(as.vector(energy),as.vector(energy.randmio),as.vector(energy.DLW))
state.actual <- rep(as.character(1:numClusters),3)
grp.actual <- factor(c(rep('SC',numClusters),rep('Deg. Pres.',numClusters),rep('Space Pres.',numClusters)),levels = c('Deg. Pres.','Space Pres.','SC'),ordered = TRUE)

# check that these vectors were made correctly

print(paste('If you see TRUE', numClusters,'times below for each null model, plotting data frames have been constructed correctly'))
print('SC')
print(sapply(1:numClusters, function(K) identical(E[state == as.character(K) & as.character(grp) == 'SC'],null_energy[,K])))
print(sapply(1:numClusters, function(K) identical(E.actual[state.actual == as.character(K) & as.character(grp.actual) == 'SC'],energy[,K])))
print('Degree preserving')
print(sapply(1:numClusters, function(K) identical(E[state == as.character(K) & as.character(grp) == 'Deg. Pres.'],null_energy.randmio[,K])))
print(sapply(1:numClusters, function(K) identical(E.actual[state.actual == as.character(K) & as.character(grp.actual) == 'Deg. Pres.'],energy.randmio[,K])))
print('Space preserving')
print(sapply(1:numClusters, function(K) identical(E[state == as.character(K) & as.character(grp) == 'Space Pres.'],null_energy.DLW[,K])))
print(sapply(1:numClusters, function(K) identical(E.actual[state.actual == as.character(K) & as.character(grp.actual) == 'Space Pres.'],energy.DLW[,K])))
energy.repmat <- kronecker(matrix(1,nperms,1),energy)
energy.mio.repmat <- kronecker(matrix(1,nperms,1),energy.randmio)
energy.DLW.repmat <- kronecker(matrix(1,nperms,1),energy.DLW)

# compute p-values, i.e. probability that persistence energy of null states is < actual states for each network type

pvals.orig <- c(p.adjust(colMeans(null_energy < energy.repmat),method='bonferroni'),  # bonferroni correct over each cluster separately for each null model / actual data
  p.adjust(colMeans(null_energy.randmio < energy.mio.repmat),method='bonferroni'),
  p.adjust(colMeans(null_energy.DLW < energy.DLW.repmat),method='bonferroni'))
pvals.size <- ifelse(pvals.orig < 0.05,yes = 2.5,no = 0)    # don't show asterisks for non-sig results by making them size =0
p.labs <- ifelse(pvals.orig < 0.05,yes = '*',no = '')       # asterisk or no asterisk
pval.y <- 1.1*c(rep(max(null_energy),numClusters),rep(max(null_energy.randmio),numClusters),rep(max(null_energy.DLW),numClusters))  # height of p-val labels

p <- ggplot() + geom_boxplot(aes(x = state,y=E,color = grp,fill=grp),position = position_dodge(width = 1),outlier.shape=20,outlier.size = 0.5,lwd = 0.3) + 
  geom_point(aes(x=state.actual,y=as.numeric(E.actual),color = grp.actual),
             shape = 24,stroke = 0,size = 1,fill = '#7F2A49',position = position_dodge(width=1)) +
  #geom_text(aes(x=state.actual[!is.na(pvals)],y=as.numeric(pval.y)[!is.na(pvals)],size=pvals[!is.na(pvals)],
   #             color=grp.actual[!is.na(pvals)]),label = '*',position = position_dodge(width=1)) +
  geom_text(aes(x=state.actual,y=as.numeric(pval.y),size=pvals.size,
                color=grp.actual,label = p.labs), position = position_dodge(width=1)) +
  scale_x_discrete(limits = 1:numClusters,breaks = 1:numClusters,labels = list(clusterNames)) +
  scale_size_continuous(guide = FALSE) + 
  scale_color_manual(limits = c('Deg. Pres.','Space Pres.','SC'),breaks = c('Deg. Pres.','Space Pres.','SC'),values = c('#005C9FFF','#71AABEFF','#E2492FFF')) +
  scale_fill_manual(limits = c('Deg. Pres.','Space Pres.','SC'),breaks = c('Deg. Pres.','Space Pres.','SC'),values = c('#005C9F1A','#71AABE1A','#E2492F1A')) +
  ylab('Min. Control Energy') + xlab('') + ggtitle('State Persistence Energy') + 
  theme_classic() + theme(text = element_text(size= 8)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8))

if(numClusters == 5 | numClusters == 6){
  p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

ggsave(plot = p, filename = paste(masterdir,'analyses/control_energy/PersistenceEnergyVsSpherePermAllNullModels_FA_k',numClusters,'c',c,'T',T,'.pdf',sep =""),
  height = 2.25,width=4, units = "in",useDingbats=FALSE)