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
restDwell <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombDwellTime_k',
                         numClusters,name_root,'.mat',sep = ''))$DwellTimeMean
nbackBlockDwell <- lapply(1:numBlocks, function(B) readMat(paste(masterdir,'analyses/nbackblocks/DwellTime',
	BlockNames[B],'_k',numClusters,name_root,'.mat',sep = ''))$BlockDwellTime)

nbackBlockDwell.dprime <- vector("list",numBlocks)	# initialize results list, extract first non-intercept coefficient (reason for index of 2)
nbackBlockDwell.dprime[[1]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBlockDwell[[1]][,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackBlockDwell.dprime[[2]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBlockDwell[[2]][,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackBlockDwell.dprime[[3]] <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackBlockDwell[[3]][,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])
restDwell.dprime <- lapply(1:numClusters, function(K) summary(lm.beta(lm(restDwell[,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# plot

pdat <- as.data.frame(cbind(t(as.data.frame(nbackBlockDwell.dprime[[1]], row.names = c('ZbackB','ZbackP'))),
	t(as.data.frame(nbackBlockDwell.dprime[[2]],row.names = c('ObackB','ObackP'))),
	t(as.data.frame(nbackBlockDwell.dprime[[3]],row.names = c('TbackB','TbackP'))),
	t(as.data.frame(restDwell.dprime,row.names = c('RestB','RestP')))))
# pdat$ZbackP <- p.adjust(pdat$ZbackP,method = 'bonf') < 0.05
# pdat$ObackP <- p.adjust(pdat$ObackP,method = 'bonf') < 0.05
# pdat$TbackP <- p.adjust(pdat$TbackP,method = 'bonf') < 0.05
# pdat$ZbackP <- ifelse(pdat$ZbackP,'*','')
# pdat$ObackP <- ifelse(pdat$ObackP,'*','')
# pdat$TbackP <- ifelse(pdat$TbackP,'*','')
p.list <- list.posthoc.correct(list(pdat$ZbackP,pdat$ObackP,pdat$TbackP,pdat$RestP),method = 'bonf')
pdat$ZbackP <- ifelse(p.list[[1]] < 0.05,'*','')
pdat$ObackP <- ifelse(p.list[[2]] < 0.05,'*','')
pdat$TbackP <- ifelse(p.list[[3]] < 0.05,'*','')
pdat$RestP <- ifelse(p.list[[4]] < 0.05,'*','')

rownames(pdat) <- clusterNames
pdat <- data.frame(B = c(pdat$RestB,pdat$ZbackB,pdat$ObackB,pdat$TbackB),
	p=c(pdat$RestP,pdat$ZbackP,pdat$ObackP,pdat$TbackP))
pdat$Scan = rep(c('00Rest','0-back','1-back','2-back'),each=numClusters)
pdat$State = rep(clusterNames,length(p.list))

BlockColors <- c("#FBB4AE", "#B3CDE3", "#CCEBC5") # brewer.pal(3,'Pastel1')
PlotColors <- c(RNcolors[1],BlockColors)
p <- ggplot(data = pdat, aes(y = B, x = State, fill = Scan,label = p)) + 
	geom_bar(stat = "identity", position = position_dodge(width = NULL),alpha = 0.8, color = 'black') +
  geom_text(position = position_dodge(width = 0.9), size = 5) + xlab("") + ylab(expression(beta["Age"])) + ggtitle('Block Dwell Time') +
  scale_fill_manual(values = PlotColors) + scale_x_discrete(limits = clusterNames,breaks = clusterNames) + 
  scale_y_continuous(limits = c(1.05*min(pdat$B),1.1*max(pdat$B)))+ 
  theme_classic() + theme(text = element_text(size = 8)) + theme(legend.key.size = unit(0.5,'line')) + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.95)) +
  theme(legend.title = element_blank(),plot.title = element_text(hjust =0.5,size = 8),legend.position = 'none') +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
if(numClusters == 5){
	p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

save(p,pdat,file=paste0(savedir,'Fig6b__RestBlockDwellTimeAge.RData'))
ggsave(plot = p,filename = paste(savedir,'RestAndnBackBlockAgeDwellTime_k',numClusters,'.pdf',sep =''),units = 'cm',height = 3,width = 6)