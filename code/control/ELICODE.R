####split data into within-person and between-person components
#This entails a couple of steps
#this gives you the average structure value for each person
DATASETNAME$AveSTRUCTURALVARIABLE <- with(DATASETNAME, ave(STRUCTUREVARIABLENAME, IDNUM, FUN=function(x) mean(x, na.rm=TRUE))) #IDNUM is your ID variable

###Get the sample mean of the between-person average of the structural variable:
describe(DATASETNAME$AveSTRUCTURALVARIABLE) #make note of the mean of this variable and plug it in below

#now we set sample mean for the structural variable

DATASETNAME$MeanSTRUCTURALVARIABLE <- #put the value found above here
  
  
#Split Structure into a between and within subjects component
  
#Subtract grand mean from the daily scores IPC "grand mean centered"
DATASETNAME$STRUCTUREc <- DATASETNAME$STRUCTURALVARIABLENAME - DATASETNAME$MeanSTRUCTURALVARIABLE

#Now we create a between-subjects mean from this grand mean centred variable
DATASETNAME$STRUCTUREbw <- with(DATASETNAME, ave(STRUCTUREc, ID, FUN = mean))

#And then a within-subjects mean using this between-subject mean
DATASETNAME$STRUCTUREwn <- DATASETNAME$STRUCTUREc - DATASETNAME$STRUCTUREbw

####so your two important variables are STRUCTUREbw and STRUCTUREwn

####Run a multilevel model
library(nlme)
options(scipen = 999)

# Predicting function

ctrl <- lmeControl(opt = 'optim')
MODELNAME <- lme(fixed=FUNCTIONALVARIABLE ~ STRUCTUREwn + STRUCTUREbw + Task + STRUCTUREwn*Task, 
                data=DATASETNAME,
                random=~ STRUCTUREwn | IDNUM, 
                na.action=na.omit, control=ctrl)
intervals(MODELNAME)
summary(MODELNAME)
vcov(MODELNAME)

