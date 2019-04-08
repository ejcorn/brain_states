args <- commandArgs(TRUE)
name_root <- 'ScanCLaus250Z0final'
numClusters <- 5
basedir <- '~/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states/'

library(ggplot2)
library(R.matlab)
library(RColorBrewer)
library(lm.beta)
library(reshape2)
library(viridis)
library(caret)

masterdir <- paste(basedir,'results/',name_root,'/',sep='')
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))
source(paste(basedir,'code/miscfxns/statfxns.R',sep=''))
savedir <- paste(masterdir,'development/SCTP/',sep = '')
dir.create(savedir,recursive = TRUE)

clusterNames <- readMat(paste(basedir,'results/',name_root,'/transitionprobabilities/OverallClusterCentroids_k',numClusters,name_root,'.mat',sep=''))
clusterNames <- unlist(clusterNames$clusterNames)
clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
RNcolors <- c('#005C9F','#FF8400') 

persistExclude <- (1:numClusters) + (numClusters*(0:(numClusters-1)));

demo <- read.csv(paste(basedir,'data/Demographics',name_root,'.csv',sep =""))

thrsh <- 1.2		# most stringent
SC <- readMat(paste(masterdir,"structrans/structrans_k",numClusters,"thresh",thrsh,name_root,".mat",sep=""))$interStateSC

restTP <- readMat(paste(masterdir,'transitionprobabilities/RestCombTransitionProbabilities_k',
	numClusters,name_root,'.mat',sep = ''))$transitionProbability
nbackTP <- readMat(paste(masterdir,'transitionprobabilities/nBackCombTransitionProbabilities_k',
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

data <- cbind(data.frame(nbackSCTPCor.r=nbackSCTPCor.r),demo)
# https://quantdev.ssri.psu.edu/tutorials/cross-validation-tutorial
# http://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/#k-fold-cross-validation
nreps <- 100
k <- 10
data_ctrl <- trainControl(method = "repeatedcv", number = k,repeats=nreps)
model_caret <- train(nbackSCTPCor.r ~ age_in_yrs + BrainSegVol + handedness + nbackRelMeanRMSMotion + Sex,   # model to fit
                     data = data,                        
                     trControl = data_ctrl,              # folds
                     method = "lm")                      # specifying regression model

full.RMSE <- summary(mdl)$sigma   # extract error from model fit on full sample
full.Rsq <- summary(mdl)$r.squared # extract R^2 from model fit on full sample
kfold.RMSE <- model_caret$resample$RMSE   # extract distribution of errors from k-fold CV
kfold.Rsq <- model_caret$resample$Rsquared   # extract distribution of R^2 from k-fold CV

# test if full RMSE differs from out-of-sample RMSE: two-tailed test asks if model is over or underfitting

print(paste('Age-r(SC,TP) n-back, p-val for RMSE in full sample vs. distribution of',nreps,'repeats of k =',k,'folds:',pval.2tail.np(test.val = full.RMSE,dist = kfold.RMSE)))
print(paste('Age-r(SC,TP) n-back, p-val for R-squared in full sample vs. distribution of',nreps,'repeats of k =',k,'folds:',pval.2tail.np(test.val = full.Rsq,dist = kfold.Rsq)))

# now for symmetry:
nreps <- 100 # set cv repetitions and number of folds
k <- 10
scanlab <- c('RestComb','nBackComb')

for(i in 1){   # only look at rest which is presented in Fig 7g
  
  # symmetry
  symmvars <- readMat(paste(basedir,'results/',name_root,'/transitionprobabilities/symmetry/RvNSymmetryScorev2_k',
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
  print(summary(age.symm))
  print(paste(scanlab[i],', ','Cohen\'s f^2 for age: ',cohens.f2(age.symm,'age_in_yrs'),sep=''))
  data <- cbind(data.frame(symm=symm),data.frame(hm),demo) 
  data_ctrl <- trainControl(method = "repeatedcv", number = k,repeats=nreps)
  model_caret <- train(symm ~ age_in_yrs + BrainSegVol + handedness + hm + Sex,   # model to fit
                       data = data,                        
                       trControl = data_ctrl,              # folds
                       method = "lm")                      # specifying regression model
  
  full.RMSE <- summary(age.symm)$sigma   # extract error from model fit on full sample
  full.Rsq <- summary(age.symm)$r.squared # extract R^2 from model fit on full sample
  kfold.RMSE <- model_caret$resample$RMSE   # extract distribution of errors from k-fold CV
  kfold.Rsq <- model_caret$resample$Rsquared   # extract distribution of R^2 from k-fold CV
  
  # test if full RMSE differs from out-of-sample RMSE: two-tailed test asks if model is over or underfitting
  
  print(paste('Age-asymm',scanlab[i],', p-val for RMSE in full sample vs. distribution of',nreps,'repeats of k =',k,'folds:',pval.2tail.np(test.val = full.RMSE,dist = kfold.RMSE)))
  print(paste('Age-asymm',scanlab[i],', p-val for R-squared in full sample vs. distribution of',nreps,'repeats of k =',k,'folds:',pval.2tail.np(test.val = full.Rsq,dist = kfold.Rsq)))
  
}
