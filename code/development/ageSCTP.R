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
savedir <- paste(masterdir,'analyses/development/SCTP/',sep = '')
dir.create(savedir,recursive = TRUE)

clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
RNcolors <- c('#005C9F','#FF8400') 

persistExclude <- (1:numClusters) + (numClusters*(0:(numClusters-1)));

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep =""))

thrsh <- 1.2		# most stringent
SC <- readMat(paste(masterdir,"analyses/structrans/structrans_k",numClusters,"thresh",thrsh,name_root,".mat",sep=""))$interStateSC

restTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/RestCombTransitionProbabilities_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability
nbackTP <- readMat(paste(masterdir,'analyses/transitionprobabilities/nBackCombTransitionProbabilities_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability

nobs <- nrow(demo)

SC <- SC[,-persistExclude]	# only look at off-diagonal elements
restTP <- restTP[,-persistExclude]
nbackTP <- nbackTP[,-persistExclude]

restSCTPCor.r <- sapply(1:nobs, function(N) cor.test(restTP[N,],SC[N,])$estimate)
nbackSCTPCor.r <- sapply(1:nobs, function(N) cor.test(nbackTP[N,],SC[N,])$estimate)

mdl <- lm.beta(lm(nbackSCTPCor.r ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, data = demo))
print(summary(mdl))
print(paste('Cohen\'s f^2 for age:',cohens.f2(mdl,'age_in_yrs')))
# construct data frame with partial residuals WRT age and age
df <- data.frame(y = residuals(mdl) + summary(mdl)$coef['age_in_yrs','Estimate']*demo$age_in_yrs, x = demo$age_in_yrs)
names(df) <- c('r(SC,TP)','Age')
pval <- signif(summary(mdl)$coef['age_in_yrs','Pr(>|t|)']*2,2) # bonferroni correct over rest and n-back
p <- ggplot(df,aes(x=Age,y=`r(SC,TP)`)) + geom_point(color = 'grey60',stroke = 0, size = 1, alpha = 0.6) + ggtitle('n-back') + 
	geom_smooth(fill = RNcolors[2], color = RNcolors[2], method = 'lm') + 
	geom_text(aes(x = 0.9*max(Age), y = 0.9*max(`r(SC,TP)`),label = paste('p =',pval)), size = 2.5) +
	theme_classic() + theme(axis.text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5,size=8)) + 
	theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
p
ggsave(plot = p, filename = paste(savedir,'n-back','SCTPVsAge_k',numClusters,'.pdf',sep =''),units = 'in', width = 2,height = 2)

mdl <- lm.beta(lm(restSCTPCor.r ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex, data = demo))
print(summary(mdl))
# construct data frame with partial residuals WRT age and age
df <- data.frame(y = residuals(mdl) + summary(mdl)$coef['age_in_yrs','Estimate']*demo$age_in_yrs, x = demo$age_in_yrs)
names(df) <- c('r(SC,TP)','Age')
pval <- signif(summary(mdl)$coef['age_in_yrs','Pr(>|t|)']*2,2)	# bonferroni correct over rest and n-back
p <- ggplot(df,aes(x=Age,y=`r(SC,TP)`)) + geom_point(color = 'grey60',stroke = 0, size = 1, alpha = 0.6) + ggtitle('Rest') + 
	geom_smooth(fill = RNcolors[2], color = RNcolors[2], method = 'lm') + 
	geom_text(aes(x = 0.9*max(Age), y = 0.9*max(`r(SC,TP)`),label = paste('p =',pval)), size = 2.5) +
	theme_classic() + theme(axis.text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5,size=8)) + 
	theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
p
ggsave(plot = p, filename = paste(savedir,'Rest','SCTPVsAge_k',numClusters,'.pdf',sep =''),units = 'in', width = 2,height = 2)