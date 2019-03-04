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

clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
RNcolors <- c('#005C9F','#FF8400')  

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep =""))

savedir <- paste(masterdir,'analyses/development/blockdur/',sep = '')
if(!dir.exists(savedir)){
	dir.create(savedir,recursive = TRUE)
}

restDwell <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombDwellTime_k',
                         numClusters,name_root,'.mat',sep = ''))$DwellTime
nbackDwell <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombDwellTime_k',
                        numClusters,name_root,'.mat',sep = ''))$DwellTime

restDwell.age <- lapply(1:numClusters, function(K) summary(lm.beta(lm(restDwell[,K] ~ age_in_yrs + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])
nbackDwell.age <- lapply(1:numClusters, function(K) summary(lm.beta(lm(nbackDwell[,K] ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, 
	data = demo)))$coefficients[2,c('Standardized','Pr(>|t|)')])

# permutation testing
nperms <- 500
age_orig <- demo$age_in_yrs
restDwell.age.Null <- lapply(1:nperms, function(n) lapply(1:numClusters, function(K) summary(lm.beta(lm(restDwell[,K] ~ age_orig + BrainSegVol + handedness + restRelMeanRMSMotion + Sex, 
  data = demo[sample(nrow(demo),replace=F),])))$coefficients[2,c('Standardized','Pr(>|t|)')]))
restBeta.null <- lapply(1:numClusters, function(K) sapply(1:nperms, function(n) restDwell.age.Null[[n]][[K]]['Standardized']))
restBeta.z <- lapply(1:numClusters, function(K) (restDwell.age[[K]]['Standardized'] - mean(restBeta.null[[K]])) /sd(restBeta.null[[K]]))
restBeta.p <- lapply(1:numClusters, function(K) restDwell.age[[K]]['Standardized'] > restBeta.null[[K]])

pdat <- as.data.frame(cbind(t(as.data.frame(restDwell.age, row.names = c('RestB','RestP'))),t(as.data.frame(nbackDwell.age,row.names = c('nBackB','nBackP')))))
pdat$RestP <- p.adjust(pdat$RestP,method = 'bonf') < 0.05
pdat$nBackP <- p.adjust(pdat$nBackP,method = 'bonf') < 0.05
pdat$RestP <- ifelse(pdat$RestP,'*','')
pdat$nBackP <- ifelse(pdat$nBackP,'*','')
colnames(pdat) <- rep(c('B','p'),2)
pdat <- as.data.frame(rbind(pdat[1:numClusters,1:2],pdat[1:numClusters,3:4]))
pdat$Scan = c(rep('Rest',numClusters),rep('n-back',numClusters))
pdat$State = rep(clusterNames,2)

p <- ggplot(data = pdat, aes(y = B, x = State, fill = Scan,label = p)) + 
	geom_bar(stat = "identity", position = position_dodge(width = NULL), color = 'black') +
  geom_text(position = position_dodge(width = 0.9), size = 5) + xlab("") + ylab(expression(beta["age"])) + ggtitle('Dwell Time') +
  scale_y_continuous(limits = c(1.05*min(pdat$B),1.1*max(pdat$B))) + 
  scale_fill_manual(limits = c('Rest','n-back'), values = RNcolors) + scale_x_discrete(limits = clusterNames,breaks = clusterNames) + 
  theme_classic() + theme(text = element_text(size = 8)) + theme(legend.position = 'none') + #theme(legend.key.size = unit(0.5,'line')) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.95)) +
  theme(legend.title = element_blank(),plot.title = element_text(hjust =0.5,size = 8)) + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

if(numClusters == 5){
	p <- p + theme(axis.text.x = element_text(color = clusterColors))
}

ggsave(plot = p,filename = paste(savedir,'RvNAgeDwellTime_k',numClusters,'.pdf',sep =''),units = 'cm',height = 3,width = 5)
