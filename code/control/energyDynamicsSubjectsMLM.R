name_root= 'ScanCLaus250Z0troubleshoot'
basedir = '~/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states/'
numClusters = 5

masterdir <- paste(basedir,'results/',name_root,'/',sep='')
savedir <- paste(masterdir,'control_energy/',sep='')

library(ggplot2)
library(R.matlab)
library(RColorBrewer)

source(paste(basedir,'code/plottingfxns/GeomSplitViolin.R',sep=''))
source(paste(basedir,'code/plottingfxns/plottingfxns.R',sep=''))

scanlab <- c("RestComb","nBackComb")
onDiag <- (1:numClusters) + (numClusters*(0:(numClusters-1)))
scanttl <- c('Rest','n-back')
RNcolors <- c('#005C9F','#FF8400')  

# compare persistence energy to persistence probabilities

matResults <- readMat(paste(savedir,'SubjectPersistenceEnergy_k',numClusters,'.mat',sep = ''))
subjectPersistenceEnergy <- matResults$subjectPersistenceEnergy

rtp <- readMat(paste(masterdir,"transitionprobabilities/",
                     scanlab[1],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability
ntp <- readMat(paste(masterdir,"transitionprobabilities/",
                     scanlab[2],'TransitionProbabilities_k',numClusters,name_root,".mat",sep = ""))$transitionProbability
rpp <- rtp[,onDiag]
npp <- ntp[,onDiag]

rpp <- readMat(paste(masterdir,"transitionprobabilities/",
                     scanlab[1],'DwellTime_k',numClusters,name_root,".mat",sep = ""))$DwellTimeMean
npp <- readMat(paste(masterdir,"transitionprobabilities/",
                     scanlab[2],'DwellTime_k',numClusters,name_root,".mat",sep = ""))$DwellTimeMean


# make matrix to assign each element of persist energy matrix/persist prob matrix to a subject
ids <- do.call('rbind',lapply(1:nrow(rpp), function(i) rep(i,numClusters)))
task1 <- do.call('rbind',lapply(1:nrow(rpp), function(i) rep(0,numClusters)))
task2 <- do.call('rbind',lapply(1:nrow(rpp), function(i) rep(1,numClusters)))

# vectorize all matrices into one dataframe
df <- data.frame(IDNUM=as.vector(ids),
                 PersistEnergy=as.vector(subjectPersistenceEnergy),
                 PersistProbs=as.vector(rpp),
                 Task=as.vector(task1))
df <- rbind(df,data.frame(IDNUM=as.vector(ids),
                          PersistEnergy=as.vector(subjectPersistenceEnergy),
                          PersistProbs=as.vector(npp),
                          Task=as.vector(task2)))

####split data into within-person and between-person components
#This entails a couple of steps
#this gives you the average structure value for each person
df$MeanPersistEnergy <- with(df, ave(PersistEnergy, IDNUM, FUN=function(x) mean(x, na.rm=TRUE))) #IDNUM is your ID variable

  #Split Structure into a between and within subjects component
  
  #Subtract grand mean from the daily scores IPC "grand mean centered"
  df$PersistEnergyCentered <- df$PersistEnergy - mean(df$MeanPersistEnergy)

#Now we create a between-subjects mean from this grand mean centred variable
df$PersistEnergyBetween <- with(df, ave(PersistEnergyCentered, IDNUM, FUN = mean))

#And then a within-subjects mean using this between-subject mean
df$PersistEnergyWithin <- df$PersistEnergyCentered - df$PersistEnergyBetween

####so your two important variables are PersistEnergyBetween and PersistEnergyWithin

####Run a multilevel model
library(nlme)
options(scipen = 999)

# Predicting function

ctrl <- lmeControl(opt = 'optim')
m <- lme(fixed=PersistProbs ~ PersistEnergyWithin + PersistEnergyBetween + Task + Task*PersistEnergyWithin, 
                 data=df,
                 random=~ PersistEnergyWithin + Task*PersistEnergyWithin | IDNUM, 
                 na.action=na.omit, control=ctrl)
intervals(m)
summary(m)
vcov(m)
