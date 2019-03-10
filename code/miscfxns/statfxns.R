pval.2tail.np <- function(test.val,dist){
	# test.val: individual value being compared to distribution
	# dist: vector,distribution of values under some null model, or otherwise
	# sig.fig: number of significant figures
	# compute 2-tailed p-value for test value occurring in distribution
	dist <- as.numeric(dist)
	pval.2tail <- 2*min(mean(test.val > dist),mean(test.val < dist))
	return(pval.2tail)
}


pval.label.np <- function(pvals,n,sig.fig=2){
	# pval: p-value obtained from non-parametric test
	# n: number of samples in distribution used to obtain p-value
	# sig.fig: number of significant figures to report
	# don't ever say p = 0; instead replace with p < 1/n
	
	# make p-value label, rounded to sig figs
	pval.txt <- paste('p =',signif(pvals,sig.fig))
	# replace p = 0 with p < 1/n
	pval.txt[pvals == 0] <- paste('p <',signif(1/n,sig.fig))
	return(pval.txt)
}