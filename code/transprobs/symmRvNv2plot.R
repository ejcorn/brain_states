# plot symmetry graphs

args <- commandArgs(TRUE)
name_root <- args[1]
numClusters <- as.numeric(args[2])
basedir <- args[3]

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

symmvars <- readMat(paste(basedir,'results/',name_root,'/analyses/transitionprobabilities/symmetry/RvNSymmetryScorev2_k',numClusters,'.mat',sep = ''))

scanlab = c('RestComb','nBackComb')
RNcolors <- c('#005C9F','#FF8400') 

#pvals <- cbind(symmvars$restSymmetryScore - symmvars$restNullSymmetryScore,symmvars$nBackSymmetryScore-symmvars$nBackNullSymmetryScore) # normalize to null
pvals <- cbind(symmvars$restSymmetryScore,symmvars$nBackSymmetryScore)	# currently leaving normalization out to make score scale more interpretable but result is same
grps <- c(rep('Rest',nrow(pvals)),rep('n-back',nrow(pvals)))

p <- ggplot() + geom_violin(aes(x = grps,y = as.vector(pvals), fill = grps)) + theme_classic() +
  scale_x_discrete(limits = c('Rest','n-back')) +
  scale_fill_manual(limits = c('Rest','n-back'), values = c(RNcolors)) +
  xlab('') + ylab('Asymmetry Score') + theme(legend.title = element_blank()) + theme(text = element_text(size = 8)) +
  theme(legend.position = 'none')
p

test.output <- t.test(symmvars$nBackSymmetryScore,symmvars$restSymmetryScore,paired = TRUE)
print(test.output)
if(test.output$p.value < 10^-15){
	p <- p + annotate("text", x = 'n-back', y = 0.5,label = "**",color = 'red')
}
ggsave(plot = p,filename = paste(basedir,'results/',
                                 name_root,'/analyses/transitionprobabilities/symmetry/RvNSymmetryScorev2_k',
                                 numClusters,name_root,'.pdf',sep = ''), units = 'in', height = 2, width = 2)