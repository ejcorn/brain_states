args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))

scanlabs <- c('Rest','nBack')  # for files
scanttls <- c('Rest','n-back') # for plot labels
for(scan in c(1:2)){
  scanlab <- scanlabs[scan]
  scanttl <- scanttls[scan]

  pncDwell <- readMat(paste(masterdir,'analyses/transitionprobabilities/',scanlab,'CombDwellTime_k',numClusters,name_root,'.mat',sep = ''))$DwellTimeMean
  hcpDwell <- readMat(paste(masterdir,'analyses/hcpLR/',scanlab,'CombHCPDwellTime_k',numClusters,name_root,'.mat',sep = ''))$HCPDwellTimeMean
  clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
  clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])
  clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
  PHcolors <- c('#005C9F','#33C3A6')  

  grps <- rbind(matrix('PNC',nrow = nrow(pncDwell),ncol = ncol(pncDwell)),matrix('HCP',nrow = nrow(hcpDwell),ncol = ncol(hcpDwell)))
  statesPNC <- sapply(1:numClusters, function(K) rep(as.character(K),nrow(pncDwell)))
  statesHCP <- sapply(1:numClusters, function(K) rep(as.character(K),nrow(hcpDwell)))
  df <- data.frame(grps =as.vector(grps),states=as.vector(rbind(statesPNC,statesHCP)),dwelltimes=as.vector(rbind(pncDwell,hcpDwell)))
  p <- ggplot(data=df) + geom_split_violin(aes(x = states ,y = dwelltimes,fill = grps)) + theme_classic() +
    scale_fill_manual(limits = c('PNC','HCP'), values = PHcolors) + 
    ylab("Dwell Time (seconds)") + xlab("") +theme(text = element_text(size = 8)) +
    theme(legend.title = element_blank()) + 
    scale_x_discrete(limits = 1:numClusters, breaks=1:numClusters, labels = list(clusterNames)) +
    theme(legend.key.size = unit(0.5,'line'))

  diffs.full <- lapply(1:numClusters, function(i) t.test(pncDwell[,i],hcpDwell[,i],paired=FALSE))
  diffs <- sapply(diffs.full, function(x) x$p.value)
  diffs <- p.adjust(diffs,method = "bonf")
  for(K in 1:numClusters){
    if(diffs[K] < 10^-15){
      p <- p + annotate("text", x = K, y = 1.1*max(rbind(pncDwell,hcpDwell)),label = "**",color = 'red')
    } else if(diffs[K] < 10^-4){
    	p <- p + annotate("text", x = K, y = 1.1*max(rbind(pncDwell,hcpDwell)),label = "*",color = 'red')
    }
  }

  if(numClusters == 5){
    p <- p + theme(axis.text.x = element_text(size=8,colour = clusterColors))
  }

  ggsave(plot = p, filename = paste(masterdir,'analyses/hcpLR/PNCvsHCP',scanlab,'DwellTime_k',numClusters,'.pdf',sep =""),height = 2,width = (numClusters-1) + 0.25, units = "in")

}