args <- commandArgs(TRUE)
name_root <- args[1]
scan <- args[2]
numClusters <- as.numeric(args[3])
basedir <- args[4]

library(R.matlab)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

thrsh.rng <- signif(seq(0,1.5,0.1),2)
masterdir <- paste(basedir,'results/',name_root,'/',sep='')
persistExclude <- (1:numClusters) + (numClusters*(0:(numClusters-1)));

lncol <- "#CE2B37"
ptcol <- '#005C9F'
lncol.null <- 'grey50'
ptcol.null <- 'black'

if(scan == 'R'){
  scanlab = 'Rest'
} else if(scan == 'N'){
  scanlab = 'nBack'
} else if(scan == 'C'){
  scanlab = c('RestComb','nBackComb')
  scanttl = c('Rest','n-back')
}

##########################
### Raw SC vs BCT Null ###
##########################

overlap <- lapply(thrsh.rng, function(T) readMat(paste(masterdir,"analyses/structrans/",
            "structrans_k",numClusters,"thresh",T,name_root,".mat",sep = ""))$overlap)

for(S in 1:length(scanlab)){

    TP.full <- readMat(paste(masterdir,"analyses/transitionprobabilities/"
        ,scanlab[S],'TransitionProbabilities_k',numClusters,name_root,".mat",sep=""))$transitionProbability
    nobs <- nrow(TP.full)
    TP.full <- TP.full[,-persistExclude] # exclude diagonal/persistence probs
    TP <- colMeans(TP.full)
    
    SC.bythrsh <- list() # store group average SC
    SC.bythrsh.null <- list() # store group average null SC
    SCTP.subj.r <- list() # store subject level correlations
    SCTP.subj.null.r <- list() # store subject level correlations
    SCTP.r.vs.null.t <- list()   # compare subject level correlation dist. to null using t-test
    
    for(T in 1:length(thrsh.rng)){
        SC <- readMat(paste(masterdir,"analyses/structrans/",
            "structrans_k",numClusters,"thresh",thrsh.rng[T],name_root,".mat",sep = ""))$interStateSC
        SC[is.na(SC)] <- 0  #NA comes from 0 connections so make this 0
        SC <- SC[,-persistExclude] # exclude diagonal/persistence probs
        SC.bythrsh[[T]] <- colMeans(SC)

        SC.null <- readMat(paste(masterdir,"analyses/structrans/BCTnull/",
            "BCTNULLstructrans_k",numClusters,"thresh",thrsh.rng[T],name_root,".mat",sep = ""))$interStateSC
        SC.null[is.na(SC.null)] <- 0  #NA comes from 0 connections so make this 0
        SC.null <- SC.null[,-persistExclude] # exclude diagonal/persistence probs
        SC.bythrsh.null[[T]] <- colMeans(SC.null)

        # compare subject level correlations between null and actual
        r <- sapply(1:nobs, function(N) cor(SC[N,],TP.full[N,]))
        r.null <- sapply(1:nobs, function(N) cor(SC.null[N,],TP.full[N,]))
        SCTP.r.vs.null.t[[T]] <- t.test(r,r.null,paired=TRUE)
        SCTP.subj.r[[T]] <- r
        SCTP.subj.null.r[[T]] <- r.null

    }

    sig.vs.null <- sapply(SCTP.r.vs.null.t, function(x) x$estimate > 0) & sapply(SCTP.r.vs.null.t, function(x) x$p.value < 10^-15)
    
    SCTP.bythrsh <- sapply(SC.bythrsh, function(x) cor(x,TP))   # compute full sample group average SC-TP correlation
    SCTP.bythrsh.null <- sapply(SC.bythrsh.null, function(x) cor(x,TP))   # compute null sample group average SC-TP correlation
    SCTP.bythrsh <- sapply(SCTP.subj.r, function(x) mean(x))   # compute average subject level SC-TP correlation
    SCTP.bythrsh.null <- sapply(SCTP.subj.null.r, function(x) mean(x))   # compute null average subject level SC-TP correlation

    pv <- as.character(ifelse(sig.vs.null,yes='*',no=''))
    df <- data.frame(SCTP=SCTP.bythrsh,SCTP.null=SCTP.bythrsh.null,thrsh=thrsh.rng,pv=pv)#,
    save(df,file=paste(masterdir,'analyses/structrans/FigS11a-b__ThreshRange',scanlab[S],'TPSC_k',numClusters,name_root,'.RData',sep=''))

    print(scanttl[S])
    p <- ggplot(data=df) + geom_line(aes(x = thrsh,y = SCTP),color = lncol,size = 0.5) +   
        geom_point(aes(x = thrsh,y = SCTP),color = ptcol,size = 1) +
        geom_line(aes(x=thrsh,y=SCTP.null),color = lncol.null,size = 0.5) +
        geom_point(aes(x=thrsh,y=SCTP.null),color = ptcol.null, size=1) +
        geom_text(aes(x = thrsh,label = pv,y=0.3),color='red') +
        xlab('Activity Threshold (z)') + ylab('r(SC,TP)') + ggtitle(scanttl[S]) + 
        scale_x_continuous(breaks=seq(-1.5,1.5,0.3)) + scale_y_continuous(limits=c(-0.2,0.3)) +
        guides(size='none') + 
        theme_classic() + theme(text = element_text(size = 8)) + 
        theme(plot.title = element_text(size=8,hjust=0.5,face = "bold"),legend.title = element_blank(),
            legend.text = element_text(size = 8),legend.key.size = unit(0.5,'line')) +
        #theme(legend.position = 'none') + 
        theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

    ggsave(p, filename = paste(masterdir,'analyses/structrans/ThreshRange',scanlab[S],"TPSC_k",numClusters,name_root,".pdf",sep=""),units = "in", 
        width = 2.7,height = 1.5,useDingbats=FALSE)

}