library(gsubfn)

pval.2tail.np <- function(test.val,dist){
	# test.val: individual value being compared to distribution
	# dist: vector,distribution of values under some null model, or otherwise
	# sig.fig: number of significant figures
	# compute 2-tailed p-value for test value occurring in distribution
	dist <- as.numeric(dist)
	pval.2tail <- 2*min(mean(test.val >= dist),mean(test.val <= dist))
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

list.posthoc.correct <- function(X,method){
  # unlist a list, posthoc correct over all values according to "method"
  # relist the list in the same structure and return
  return(relist(flesh=p.adjust(unlist(X),method=method),skeleton=X))
}

cohens.f2 <- function(full.model,xname){
  # from Selya et al. 2012
  # f2 = (R2.ab - R2.a) / (1 - R2.ab)
  y <- full.model$model[[1]]
  covariates.ind <- which(names(full.model$model) != xname)[-1] # remove y variable index
  x.c <- full.model$model[covariates.ind]
  R2.ab<-summary(full.model)$r.sq
  R2.a<-summary(lm(y~. , data= x.c))$r.sq
  f.2 <- (R2.ab - R2.a) / (1 - R2.ab)
  return(f.2)
}