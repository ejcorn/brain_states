# analyze TP matrix properties over development for rest and n-back

args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

library(ggplot2)
library(R.matlab)
library(lm.beta)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')

print('1')

symmdir <- paste(basedir,'results/',name_root,'/analyses/development/symmetry/',sep = '')
if(!dir.exists(symmdir)){
	dir.create(symmdir,recursive = TRUE)
}

print('2')

unifdir <- paste(basedir,'results/',name_root,'/analyses/development/randnull/',sep = '')
if(!dir.exists(unifdir)){
	dir.create(unifdir,recursive = TRUE)
}

print('2')

scanlab <- c("RestComb","nBackComb")
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep =""))

clusterNames <- readMat(paste(basedir,'results/',name_root,'/clusterAssignments/k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterAssignments[[1]][[5]])

transLabels <- as.vector(sapply(1:numClusters, function(i) sapply(1:numClusters,function(j) paste(clusterNames[i],'to',clusterNames[j]))))	#label transitions
persistExclude <- (1:numClusters) + (numClusters*(0:(numClusters-1)))
transLabels <- transLabels[-persistExclude]
print(transLabels)
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")

print('defined paths, start analysis and plotting')

for(i in 1:length(scanlab)){

	# symmetry
	symmvars <- readMat(paste(basedir,'results/',name_root,'/analyses/transitionprobabilities/symmetry/RvNSymmetryScorev2_k',
		numClusters,'.mat',sep = ''))
	if(i == 1){
		symm <- symmvars$restSymmetryScore
		hm <- demo$restRelMeanRMSMotion
	}
	if(i == 2){
		symm <- symmvars$nBackSymmetryScore
		hm <- demo$nbackRelMeanRMSMotion
	}

	age.symm <- lm.beta(lm(symm ~ age_in_yrs + BrainSegVol + handedness + hm + Sex, data = demo))
	pval <- signif(summary(age.symm)$coefficients['age_in_yrs','Pr(>|t|)'],2)
	# get partial residuals
	Asymmetry <- residuals(age.symm) + summary(age.symm)$coefficients['age_in_yrs','Estimate']*demo$age_in_yrs
	Age <- demo$age_in_yrs	
	
	p <- ggplot() + geom_point(aes(x = Age, y = Asymmetry),color = 'grey60',stroke = 0, size = 1, alpha = 0.6) + 
	geom_smooth(aes(x=Age, y = Asymmetry),fill = RNcolors[i], color = RNcolors[i], method = 'lm') + ggtitle(paste(scanttl[i])) + 
	geom_text(aes(x = 0.9*max(Age), y = 0.9*max(Asymmetry),label = paste('p =',pval)), size = 2.5) +
	theme_classic() + theme(axis.text = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5,size=8)) + 
	theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

	ggsave(plot = p, filename = paste(symmdir, scanlab[i],'SymmVsAge_k',numClusters,'.pdf',sep =''),units = 'in', width = 2,height = 2)
	# non-uniformity


}