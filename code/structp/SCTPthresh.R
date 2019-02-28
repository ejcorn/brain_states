args <- commandArgs(TRUE)
name_root <- args[1]
scan <- args[2]
numClusters <- as.numeric(args[3])

library(R.matlab)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

thrsh.rng <- signif(seq(0,1.5,0.1),2)
masterdir <- paste("/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_",name_root,"/",sep="")
persistExclude <- (1:numClusters) + (numClusters*(0:(numClusters-1)));

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

    TP.full <- readMat(paste(masterdir,"analyses/transitionprobabilities/"
        ,scanlab[S],'TransitionProbabilities_k',numClusters,name_root,".mat",sep=""))$transitionProbability
    nobs <- nrow(TP.full)
    TP.full <- TP.full[,-persistExclude] # exclude diagonal/persistence probs
    TP <- colMeans(TP.full)
    SCTP.bythrsh <- list()
    for(T in 1:length(thrsh.rng)){
        SC <- readMat(paste(masterdir,"analyses/structrans/",
            "structrans_k",numClusters,"thresh",thrsh.rng[T],name_root,".mat",sep = ""))$interStateSC
        SC[is.na(SC)] <- 0  #NA comes from 0 connections so make this 0
        SC <- colMeans(SC)[-persistExclude] # exclude diagonal/persistence probs
        SCTP.bythrsh[[T]] <- cor.test(SC,TP)
    }

    SCTP.r.vs.null <- list() #matrix(nrow=length(thrsh.rng))    # compare subject level correlation dist. to null
    for(T in 1:length(thrsh.rng)){
        SC <- readMat(paste(masterdir,"analyses/structrans/",
            "structrans_k",numClusters,"thresh",thrsh.rng[T],name_root,".mat",sep = ""))$interStateSC
        SC[is.na(SC)] <- 0  #NA comes from 0 connections so make this 0
        SC <- SC[,-persistExclude]
        
        SC.null <- readMat(paste(masterdir,"analyses/structrans/BCTnull/",
            "BCTNULLstructrans_k",numClusters,"thresh",thrsh.rng[T],name_root,".mat",sep = ""))$interStateSC
        SC.null[is.na(SC.null)] <- 0  #NA comes from 0 connections so make this 0
        SC.null <- SC.null[,-persistExclude] # exclude diagonal/persistence probs
        
        r <- sapply(1:nobs, function(N) cor(SC[N,],TP.full[N,]))
        r.null <- sapply(1:nobs, function(N) cor(SC.null[N,],TP.full[N,]))
        SCTP.r.vs.null[[T]] <- t.test(r,r.null,paired=TRUE)
    }
    sig.vs.null <- sapply(SCTP.r.vs.null, function(x) x$estimate > 0) & sapply(SCTP.r.vs.null, function(x) x$p.value < 10^-15)
    pv <- sapply(SCTP.bythrsh, function(x) x$p.value)   
    pv <- ifelse(pv<0.05,yes='p < 0.05',no='p > 0.05')  
    SCTP.bythrsh <- sapply(SCTP.bythrsh, function(x) x$estimate)
    print(scanttl[S])

    p <- ggplot() + geom_line(aes(x = thrsh.rng,y = SCTP.bythrsh),color = ptcol,size = 0.5) +   
        geom_point(aes(x = thrsh.rng,y = SCTP.bythrsh),color = '#005C9F',size = 1) +
        geom_text(aes(x = thrsh.rng[sig.vs.null],y=0.7),label = '*',color='red') +
        #scale_color_manual(limits = c('p < 0.05','p > 0.05'),values=c('#71AABE','#005C9F')) +
        xlab('Activity Threshold (z)') + ylab('r(SC,TP)') + ggtitle(scanttl[S]) + 
        scale_x_continuous(breaks=seq(-1.5,1.5,0.3)) + scale_y_continuous(limits=c(-0.3,0.7)) +
        guides(size='none') + 
        theme_classic() + theme(text = element_text(size = 8)) + 
        theme(plot.title = element_text(size=8,hjust=0.5,face = "bold"),legend.title = element_blank(),
            legend.text = element_text(size = 8),legend.key.size = unit(0.5,'line')) +
        #theme(legend.position = 'none') + 
        theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

    ggsave(p, filename = paste(masterdir,'analyses/structrans/ThreshRange',scanlab[S],"TPSC_k",numClusters,name_root,".pdf",sep=""),units = "in", width = 2.7,height = 1.5)

}


###################
### BCT Null SC ###
###################

for(S in 1:length(scanlab)){

    TP <- readMat(paste(masterdir,"analyses/transitionprobabilities/"
        ,scanlab[S],'TransitionProbabilities_k',numClusters,name_root,".mat",sep=""))$transitionProbability
    TP <- colMeans(TP)
    TP <- TP[-persistExclude] # exclude diagonal/persistence probs
    SCTP.bythrsh <- list()
    for(T in 1:length(thrsh.rng)){
        SC <- readMat(paste(masterdir,"analyses/structrans/BCTnull/",
            "BCTNULLstructrans_k",numClusters,"thresh",thrsh.rng[T],name_root,".mat",sep = ""))$interStateSC
        SC[is.na(SC)] <- 0  #NA comes from 0 connections so make this 0
        SC <- colMeans(SC)[-persistExclude] # exclude diagonal/persistence probs
        SCTP.bythrsh[[T]] <- cor.test(SC,TP)
    }
    
    pv <- sapply(SCTP.bythrsh, function(x) x$p.value)   
    pv <- ifelse(pv<0.05,yes='p < 0.05',no='p > 0.05')  
    SCTP.bythrsh <- sapply(SCTP.bythrsh, function(x) x$estimate)
    print(scanttl[S])

    p <- ggplot() + geom_line(aes(x = thrsh.rng,y = SCTP.bythrsh),color = ptcol,size = 0.5) +   
        geom_point(aes(x = thrsh.rng,y = SCTP.bythrsh),color = '#005C9F',size = 1) +  
        #scale_color_manual(limits = c('p < 0.05','p > 0.05'),values=c('#71AABE','#005C9F')) +
        xlab('Activity Threshold (z)') + ylab('r(SC,TP)') + ggtitle(paste(scanttl[S],'Null')) + 
        scale_x_continuous(breaks=seq(-1.5,1.5,0.3)) + scale_y_continuous(limits=c(-0.3,0.7)) +
        guides(size='none') + 
        theme_classic() + theme(text = element_text(size = 8)) + 
        theme(plot.title = element_text(size=8,hjust=0.5,face = "bold"),legend.title = element_blank(),
            legend.text = element_text(size = 8),legend.key.size = unit(0.5,'line')) +
        #theme(legend.position = 'none') + 
        theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

    ggsave(p, filename = paste(masterdir,'analyses/structrans/BCTnull/BCTNULLThreshRange',scanlab[S],"TPSC_k",numClusters,name_root,".pdf",sep=""),units = "in", width = 2.7,height = 1.5)

}

####################
### Rick Null SC ###
####################

for(S in 1:length(scanlab)){

    TP <- readMat(paste(masterdir,"analyses/transitionprobabilities/"
        ,scanlab[S],'TransitionProbabilities_k',numClusters,name_root,".mat",sep=""))$transitionProbability
    TP <- colMeans(TP)
    TP <- TP[-persistExclude] # exclude diagonal/persistence probs
    SCTP.bythrsh <- list()
    for(T in 1:length(thrsh.rng)){
        SC <- readMat(paste(masterdir,"analyses/structrans/",
            "NULLstructrans_k",numClusters,"thresh",thrsh.rng[T],name_root,".mat",sep = ""))$interStateSC
        SC[is.na(SC)] <- 0  #NA comes from 0 connections so make this 0
        SC <- colMeans(SC)[-persistExclude] # exclude diagonal/persistence probs
        SCTP.bythrsh[[T]] <- cor.test(SC,TP)
    }
    pv <- sapply(SCTP.bythrsh, function(x) x$p.value)   
    pv <- ifelse(pv<0.05,yes='p < 0.05',no='p > 0.05')  
    SCTP.bythrsh <- sapply(SCTP.bythrsh, function(x) x$estimate)
    print(scanttl[S])

    p <- ggplot() + geom_line(aes(x = thrsh.rng,y = SCTP.bythrsh),color = ptcol,size = 0.5) +   
        geom_point(aes(x = thrsh.rng,y = SCTP.bythrsh),color = '#005C9F',size = 1) +  
        #scale_color_manual(limits = c('p < 0.05','p > 0.05'),values=c('#71AABE','#005C9F')) +
        xlab('Activity Threshold (z)') + ylab('r(SC,TP)') + ggtitle(scanttl[S]) + 
        scale_x_continuous(breaks=seq(-1.5,1.5,0.3)) + scale_y_continuous(limits=c(-0.3,0.7)) +
        guides(size='none') + 
        theme_classic() + theme(text = element_text(size = 8)) + 
        theme(plot.title = element_text(size=8,hjust=0.5,face = "bold"),legend.title = element_blank(),
            legend.text = element_text(size = 8),legend.key.size = unit(0.5,'line')) +
        #theme(legend.position = 'none') + 
        theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

    ggsave(p, filename = paste(masterdir,'analyses/structrans/NULLThreshRange',scanlab[S],"TPSC_k",numClusters,name_root,".pdf",sep=""),units = "in", width = 2.7,height = 1.5)

}