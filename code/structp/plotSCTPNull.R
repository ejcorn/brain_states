args <- commandArgs(TRUE)
name_root <- args[1]
scan <- args[2]
thrsh <- as.numeric(args[3])
numClusters <- as.numeric(args[4])
basedir <- args[5]

library(R.matlab)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')
persistExclude <- (1:numClusters) + (numClusters*(0:(numClusters-1)));

SC <- readMat(paste(masterdir,"analyses/structrans/",
	"NULLstructrans_k",numClusters,"thresh",thrsh,name_root,".mat",sep = ""))$interStateSC
SC[is.na(SC)] <- 0	#NA comes from 0 connections so make this 0
SC <- colMeans(SC)[-persistExclude] # exclude diagonal/persistence probs

lncol <- "#CE2B37"
ptcol <- "#CE2B37"

if(scan == 'R'){
  scanlab = 'Rest'
} else if(scan == 'N'){
  scanlab = 'nBack'
} else if(scan == 'C'){
  scanlab = c('RestComb','nBackComb')
  scanttl = c('Rest','n-back')
}

##############
### Raw SC ###
##############
for(S in 1:length(scanlab)){

	TP <- readMat(paste(masterdir,"analyses/transitionprobabilities/"
		,scanlab[S],'TransitionProbabilities_k',numClusters,name_root,".mat",sep=""))$transitionProbability
	TP <- colMeans(TP)
	TP <- TP[-persistExclude] # exclude diagonal/persistence probs

	SCTP.trans.r <- cor.test(SC,TP)		
	r.trans <- SCTP.trans.r$estimate
	p.trans <- SCTP.trans.r$p.value

	p <- ggplot() + geom_point(aes(x = SC,y = TP),color = ptcol,size = 1, alpha = 0.75,stroke = 0) + 
		geom_smooth(aes(x = SC,y = TP),color = lncol, fill = lncol,method = 'lm',size=0.5) +
		scale_y_continuous(limits = c(0.05,0.26),breaks= seq(0.05,0.25,length.out=5)) + scale_x_continuous() +
		annotate("text",size = 2, x = 0.7*max(SC),y = 0.08, label = paste("r = ",signif(r.trans,2))) +
		annotate("text",size = 2, x = 0.7*max(SC),y = 0.06, label = paste("p = ",signif(p.trans,2))) +
		xlab('Mean SC') + ylab('Transition Probability') + ggtitle(scanttl[S]) + 
		theme_classic() + theme(text = element_text(size = 8)) + 
		theme(plot.title = element_text(size=8,hjust=0.5,face = "bold")) +
		theme(legend.position = 'none') + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

	ggsave(p, filename = paste(masterdir,'analyses/structrans/NULL',scanlab[S],"TPSC_k",numClusters,"thresh",thrsh,name_root,".pdf",sep=""),units = "in", width = 1.5,height = 1.5)
}

#######################
### Communicability ###
#######################

SC <- readMat(paste(masterdir,"analyses/structrans/",
	"NULLstructrans_k",numClusters,"thresh",thrsh,name_root,".mat",sep = ""))$interStateComm
SC[is.na(SC)] <- 0	#NA comes from 0 connections so make this 0
SC <- colMeans(SC)[-persistExclude] # exclude diagonal/persistence probs

for(S in 1:length(scanlab)){

	TP <- readMat(paste(masterdir,"analyses/transitionprobabilities/"
		,scanlab[S],'TransitionProbabilities_k',numClusters,name_root,".mat",sep=""))$transitionProbability
	TP <- colMeans(TP)
	TP <- TP[-persistExclude] # exclude diagonal/persistence probs

	SCTP.trans.r <- cor.test(SC,TP)		
	r.trans <- SCTP.trans.r$estimate
	p.trans <- SCTP.trans.r$p.value

	p <- ggplot() + geom_point(aes(x = SC,y = TP),color = ptcol,size = 1, alpha = 0.75,stroke = 0) + 
		geom_smooth(aes(x = SC,y = TP),color = lncol, fill = lncol,method = 'lm',size=0.5) +
		scale_y_continuous(limits = c(0.05,0.26),breaks= seq(0.05,0.25,length.out=5)) + scale_x_continuous() +
		annotate("text",size = 2, x = 0.7*max(SC),y = 0.08, label = paste("r = ",signif(r.trans,2))) +
		annotate("text",size = 2, x = 0.7*max(SC),y = 0.06, label = paste("p = ",signif(p.trans,2))) +
		xlab('Communicability') + ylab('Transition Probability') + ggtitle(scanttl[S]) + 
		theme_classic() + theme(text = element_text(size = 8)) + 
		theme(plot.title = element_text(size=8,hjust=0.5,face = "bold")) +
		theme(legend.position = 'none') + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

	ggsave(p, filename = paste(masterdir,'analyses/structrans/NULL',scanlab[S],"TPComm_k",numClusters,"thresh",thrsh,name_root,".pdf",sep=""),units = "in", width = 1.5,height = 1.5)
}