args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
thrsh <- as.numeric(args[3])
basedir <- args[4]

print(paste('Threshold:',thrsh))

library(R.matlab)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

strucdir <- paste(basedir,"results/",name_root,"/analyses/structrans/",sep = "")
transdir <- paste(basedir,"results/",name_root,"/analyses/transitionprobabilities/",sep = "")
RNcolors <- c('#005C9F','#FF8400') 

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))

scanlab = c('RestComb','nBackComb')
persistExclude <- (1:numClusters) + (numClusters*(0:(numClusters-1)));

rTP <- readMat(paste(transdir,scanlab[1],'TransitionProbabilities_k',numClusters,name_root,".mat",sep=""))$transitionProbability[,-persistExclude]
nTP <- readMat(paste(transdir,scanlab[2],'TransitionProbabilities_k',numClusters,name_root,".mat",sep=""))$transitionProbability[,-persistExclude]

structure <- readMat(paste(strucdir,"structrans_k",numClusters,"thresh",as.character(thrsh),name_root,".mat",sep=""))
structure$interStateSC[is.na(structure$interStateSC)] <- 0  #NA comes from 0 connections so make this 0
SC <- structure$interStateSC[,-persistExclude]
STP <- (structure$interStateSTP*(structure$interStateSTP < Inf))[,-persistExclude]
structure$interStateComm[is.na(structure$interStateComm)] <- 0  #NA comes from 0 connections so make this 0
Comm <- structure$interStateComm[,-persistExclude]

nobs <- nrow(SC)
  
rTPSC <- sapply(1:nobs, function(i) cor(rTP[i,],SC[i,], use = "pairwise.complete.obs"))
rTPSTP <- sapply(1:nobs, function(i) cor(rTP[i,],STP[i,], use = "pairwise.complete.obs"))
rTPComm <- sapply(1:nobs, function(i) cor(rTP[i,],Comm[i,], use = "pairwise.complete.obs"))
nTPSC <- sapply(1:nobs, function(i) cor(nTP[i,],SC[i,], use = "pairwise.complete.obs"))
nTPSTP <- sapply(1:nobs, function(i) cor(nTP[i,],STP[i,], use = "pairwise.complete.obs"))
nTPComm <- sapply(1:nobs, function(i) cor(nTP[i,],Comm[i,], use = "pairwise.complete.obs"))

print('SC, rest')
print(summary(rTPSC))
print('STP, rest')
print(summary(rTPSTP))
print('Comm, rest')
print(summary(rTPComm))
print('SC, n-back')
print(summary(nTPSC))
print('STP, n-back')
print(summary(nTPSTP))
print('Comm, n-back')
print(summary(nTPComm))

rest <- cbind(rTPSC,rTPSTP,rTPComm)
nback <- cbind(nTPSC,nTPSTP,nTPComm)

grps <- rbind(matrix('Rest',nrow = nrow(rest),ncol = ncol(rest)),matrix('n-back',nrow = nrow(nback),ncol = ncol(nback)))
metrics <- sapply(1:ncol(rest), function(K) rep(as.character(K),nrow(rest)))
df.plt <- data.frame(metrics.cat=as.vector(rbind(metrics,metrics)),r.vals=as.vector(rbind(rest,nback)),grps.cat=as.vector(grps))

if(name_root == 'ScanCLaus250Z0final' & numClusters == 5){
  save(rest,nback,df.plt,file=paste(strucdir,'Fig5d__RvNSubjStrucTPCorrs_Thresh',thrsh,'_k',numClusters,name_root,'.RData',sep = ''))
} else if(name_root == 'ScanCLaus250Z0final' & numClusters == 6){
  save(rest,nback,df.plt,file=paste(strucdir,'FigS14g__RvNSubjStrucTPCorrs_Thresh',thrsh,'_k',numClusters,name_root,'.RData',sep = ''))
} else if(name_root == 'ScanCLaus250Z0cosinefinal' & numClusters == 5){
  save(rest,nback,df.plt,file=paste(strucdir,'FigS12g__RvNSubjStrucTPCorrs_Thresh',thrsh,'_k',numClusters,name_root,'.RData',sep = ''))
} else if(name_root == 'ScanCLaus125Z0final' & numClusters == 5){
  save(rest,nback,df.plt,file=paste(strucdir,'FigS15g__RvNSubjStrucTPCorrs_Thresh',thrsh,'_k',numClusters,name_root,'.RData',sep = ''))
}

p <- ggplot(df.plt) + geom_split_violin(aes(x = metrics.cat, y = r.vals,fill = grps.cat)) + theme_classic() +
  scale_fill_manual(limits = c('Rest','n-back'), values = RNcolors) + scale_y_continuous(limits = c(-0.9,0.9)) +
  ylab("r(SC,TP)") + xlab("") +theme(text = element_text(size = 8)) +
  theme(legend.title = element_blank(),legend.position = 'none',axis.text.x = element_text(color = c("#D12631","#FFDC01","#0063B9"))) + 
  scale_x_discrete(limits = 1:ncol(rest), breaks=1:ncol(rest), labels = list(c('SC','STP','Communicability'))) 

rvn <- sapply(1:ncol(rest), function(i) t.test(rest[,i],nback[,i],paired = TRUE)$p.value)
rvn.t <- lapply(1:ncol(rest), function(i) t.test(rest[,i],nback[,i],paired = TRUE))
rv0.t <- lapply(1:ncol(rest), function(i) t.test(rest[,i])) 
nv0.t <- lapply(1:ncol(nback), function(i) t.test(nback[,i])) 
rv0 <- sapply(1:ncol(rest), function(i) t.test(rest[,i])$p.value)
nv0 <- sapply(1:ncol(nback), function(i) t.test(nback[,i])$p.value)
for(K in 1:ncol(rest)){
  if(rvn[K] < 10^-15){
    p <- p + annotate("text", x = K, y = 0.9,label = "**",color = 'red')
  } else if(rvn[K] < 10^-4){
    p <- p + annotate("text", x = K, y = 0.9,label = "*",color = 'red')
  }
  if(rv0[K] < 10^-15 & nv0[K] < 10^-15){
  	p <- p + annotate("text", x = K, y = 0.8,label = "**", color = 'blue')
  } else if(rv0[K] < 10^-6 & nv0[K] < 10^-6){
    p <- p + annotate("text", x = K, y = 0.8,label = "*", color = 'blue')
  }

}

print(rvn.t)
ggsave(plot = p, filename = paste(strucdir,'RvNStrucTPCorr_k',numClusters,'thresh',thrsh,'.pdf',sep =""),height = 4,width = 4, units = "cm")