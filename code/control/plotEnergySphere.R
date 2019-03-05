args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])

name_root <- 'ScanCLaus250Z0corrfinal'
numClusters <- 5
masterdir <- paste("/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_",name_root,"/",sep="")
masterdir <- paste("~/Dropbox/Cornblath_Bassett_Projects/code/control_fc/restnbackpipeline/results/",name_root,"/",sep="")

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

scanlab <- c("RestComb","nBackComb")
persistExclude = (1:numClusters) + (numClusters*(0:(numClusters-1)))
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  
clusterNames <- c("DMN+", "DMN-", "FPN+", "VIS+", "VIS-")
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")

source('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode/GeomSplitViolin.R')
#source('~/Dropbox/Cornblath_Bassett_Projects/code/control_fc/restnbackpipeline/analysiscode/GeomSplitViolin.R')
#source('~/Dropbox/Cornblath_Bassett_Projects/code/control_fc/restnbackpipeline/analysiscode/plottingfxns.R')

savedir = paste(masterdir,'/analyses/control_energy',sep='');
dir.create(savedir,recursive = TRUE)

# compare energy of maintaining each state with sphere permuted states for actual network and all null models

energy <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k_',numClusters,'.mat',sep = ''))$Epersist
null_energy <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k_',numClusters,'.mat',sep = ''))$Epersist.Null
energy.randmio <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k_',numClusters,'.mat',sep = ''))$Epersist.mio
null_energy.randmio <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k_',numClusters,'.mat',sep = ''))$Epersist.Null.mio
energy.DLW <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k_',numClusters,'.mat',sep = ''))$Epersist.DLW
null_energy.DLW  <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k_',numClusters,'.mat',sep = ''))$Epersist.Null.DLW

nperms <- nrow(null_energy)
state <- rep(as.vector(sapply(1:numClusters, function(K) rep(as.character(K),nperms))),3)
E <- c(as.vector(null_energy),as.vector(null_energy.randmio),as.vector(null_energy.DLW))
grp <- factor(c(rep('SC',nperms*numClusters),rep('Deg. Pres.',nperms*numClusters),rep('Space Pres.',nperms*numClusters)),levels = c('Deg. Pres.','Space Pres.','SC'),ordered = TRUE)
E.actual <- c(as.vector(energy),as.vector(energy.randmio),as.vector(energy.DLW))
state.actual <- rep(as.character(1:numClusters),3)
grp.actual <- factor(c(rep('SC',numClusters),rep('Deg. Pres.',numClusters),rep('Space Pres.',numClusters)),levels = c('Deg. Pres.','Space Pres.','SC'),ordered = TRUE)

energy.repmat <- kronecker(matrix(1,nperms,1),energy)
energy.mio.repmat <- kronecker(matrix(1,nperms,1),energy.randmio)
energy.DLW.repmat <- kronecker(matrix(1,nperms,1),energy.DLW)
pvals.orig <- c(colMeans(null_energy < energy.repmat),colMeans(null_energy.randmio < energy.mio.repmat),colMeans(null_energy.DLW < energy.DLW.repmat))
pvals <- ifelse(pvals.orig < 0.05,yes = 2.5,no = 0)
p.labs <- ifelse(pvals.orig < 0.05,yes = '*',no = '')
pval.y <- 1.1*c(rep(max(null_energy),numClusters),rep(max(null_energy.randmio),numClusters),rep(max(null_energy.DLW),numClusters))
p <- ggplot() + geom_boxplot(aes(x = state,y=E,color = grp,fill=grp),position = position_dodge(width = 1),outlier.shape=20,outlier.size = 0.5,lwd = 0.3) + 
  geom_point(aes(x=state.actual,y=as.numeric(E.actual),color = grp.actual),
             shape = 24,stroke = 0,size = 1,fill = '#7F2A49',position = position_dodge(width=1)) +
  #geom_text(aes(x=state.actual[!is.na(pvals)],y=as.numeric(pval.y)[!is.na(pvals)],size=pvals[!is.na(pvals)],
   #             color=grp.actual[!is.na(pvals)]),label = '*',position = position_dodge(width=1)) +
  geom_text(aes(x=state.actual,y=as.numeric(pval.y),size=pvals,
                color=grp.actual,label = p.labs), position = position_dodge(width=1)) +
  scale_x_discrete(limits = 1:numClusters,breaks = 1:numClusters,labels = clusterNames) +
  scale_size_continuous(guide = FALSE) + 
  scale_color_manual(limits = c('Deg. Pres.','Space Pres.','SC'),breaks = c('Deg. Pres.','Space Pres.','SC'),values = c('#005C9FFF','#71AABEFF','#E2492FFF')) +
  scale_fill_manual(limits = c('Deg. Pres.','Space Pres.','SC'),breaks = c('Deg. Pres.','Space Pres.','SC'),values = c('#005C9F1A','#71AABE1A','#E2492F1A')) +
  ylab('Min. Control Energy') + xlab('') + ggtitle('State Persistence Energy') + 
  theme_classic() + theme(text = element_text(size= 8)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 8))

if(numClusters == 5){
  p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

p

ggsave(plot = p, filename = paste(masterdir,'analyses/control_energy/PersistenceEnergyVsSpherePermAllNullModels_k',numClusters,'.pdf',sep =""),height = 3,width=4, units = "in")

# compare energy of maintaining each state with sphere permuted states for actual networks only
energy <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k_',numClusters,'.mat',sep = ''))$Epersist
null_energy <- readMat(paste(masterdir,'analyses/control_energy/PersistEnergySpherePerm_k_',numClusters,'.mat',sep = ''))$Epersist.Null
nperms <- nrow(null_energy)
state <- as.vector(sapply(1:numClusters, function(K) rep(as.character(K),nperms)))
E <- as.vector(null_energy)
grp <- rep('SC',nperms*numClusters)
E.actual <- as.vector(energy)
state.actual <- rep(as.character(1:numClusters),1)
grp.actual <- rep('SC',numClusters)
p <- ggplot() + geom_boxplot(aes(x = state,y=E),color = '#005C9FFF',fill = '#005C9F1A',position = position_dodge(width = 1),outlier.shape=20,outlier.size = 0.5,lwd = 0.5) + 
  geom_point(aes(x=state.actual,y=as.numeric(E.actual),color = grp.actual),
             shape = 24,stroke = 0,size = 1,fill = '#CE2B37',position = position_dodge(width=1)) +
  scale_color_manual(limits = c('SC','Null'),breaks = c('SC','Null'),values = c('#005C9F','#71AABE')) +
  scale_x_discrete(limits = 1:numClusters,breaks = 1:numClusters,labels = clusterNames) +
  ylab('Min. Control Energy') + xlab('') + ggtitle('Min. State Maintenance Energy') + 
  theme_classic() + theme(text = element_text(size= 8)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'none')
p
energy.repmat <- kronecker(matrix(1,nperms,1),energy)
pvals <- colMeans(null_energy < energy.repmat)
col <- rep('black',numClusters)
col[pvals < 0.05] <- 'red'
for(K in 1:numClusters){
  p <- p + annotate("text", x = as.character(K), y = 1.1*max(c(energy,null_energy)),label = paste('p =',signif(pvals[K],2)),color = col[K])
}

if(numClusters == 5){
  p <- p + theme(axis.text.x = element_text(color = clusterColors))
}
p

ggsave(plot = p, filename = paste(masterdir,'analyses/control_energy/PersistenceEnergyVsSpherePermSCOnly_k',numClusters,'.pdf',sep =""),height = 2,width=3, units = "in")
