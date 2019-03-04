#!/bin/bash
set -euxo pipefail

[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

######################
### Specify inputs ###
######################

# set parcellation, distance method, z-score, scans

ZDIM=0
METHOD='correlation'
NPARC=250
SCAN='C'
LAB='corrfinal'

ROOT='Scan'$SCAN'Laus'$NPARC'Z'$ZDIM$LAB
BASEDIR='/data/tesla-data/ecornblath/brain_states/'
MATPATH=/share/apps/matlab/R2017a/bin/matlab
RPATH=/share/apps/R/R-3.2.5/bin/Rscript
MASTERDIR=$BASEDIR'results/'$ROOT
if [ ! -d "$MASTERDIR" ]; then
  mkdir -p $MASTERDIR		# recursively create general results folder and output folder
fi

cd $BASEDIR'jobs'		# change to directory containing all shell scripts

####################
### Process data ###
####################


PROC=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q -v ZDIM=$ZDIM,NPARC=$NPARC,SCAN=$SCAN,LAB=$LAB,BD=$BASEDIR,MP=$MATPATH ProcessDataBaumSample.sh)
PROC="${PROC//[!0-9]/}"
NULLSC=$(qsub -l h_vmem=30.5G,s_vmem=30G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,NPARC=$NPARC,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC nullSC.sh)
NULLSC="${NULLSC//[!0-9]/}"

##########################
### k-means clustering ###
##########################

NREPS=1
NSPLITS=10
REPK='kmeans'$ROOT
SPLITHALVES='splithalves'$ROOT
BEGINCOMMENT
for K in {2..18}
do
	for S in $(seq $NSPLITS)
	do
	qsub -N "$REPK" -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,S=$S,D=$ROOT,N=$NREPS,Z=$ZDIM,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC kmeans.sh
	done
done

#########################
### Assess clustering ###
#########################

K=5

ASSIGN0=$(qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,BD=$BASEDIR,MP=$MATPATH -hold_jid "$REPK" getassignmentspy.sh)
ASSIGN0="${ASSIGN0//[!0-9]/}"
ASSIGN=$(qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "$ASSIGN0" reorderClusters.sh)
ASSIGN="${ASSIGN//[!0-9]/}"


for K in {2..18}
do
qsub -N "zrand" -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN calczrand.sh
done
ENDCOMMENT
ASSIGN={}
qsub -N "plotzrand" -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,BD=$BASEDIR,MP=$MATPATH -hold_jid "zrand" plotzrand.sh

qsub -l h_vmem=10.5G,s_vmem=10G -q all.q,basic.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN stateRepresentation.sh

################
### Analysis ###
################

NPERMS=5000
TP='transprobs'
SYMM='symm'
K=5

######################
### State Dynamics ###
######################

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN getcentroids.sh
TP=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN getDynamics.sh)
TP="${TP//[!0-9]/}"

NULL=$(qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,N=$NPERMS,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP nullTransProbs.sh)
NULL="${NULL//[!0-9]/}"
qsub -N "$SYMM" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP symmRvNv2.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=MATPATH -hold_jid $TP persistNull.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=MATPATH -hold_jid $NULL restvsnbackNP.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP plotDur.sh

qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=MATPATH,RP=$RPATH -hold_jid $TP nbackblockDur.sh
NBACK="${NBACK//[!0-9]/}"

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP TPDistance.sh

###################
### Development ###
###################

qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,RP=$RP -hold_jid $TP ageTPDur.sh
qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,RP=$RP -hold_jid "$SYMM" ageTPprop.sh

##########################
### Structure-function ###
##########################
STRUC='struc'$ROOT
STRUCNULL='strucnull'$ROOT
BCTSTRUCNULL='BCTstrucnull'$ROOT

for THRESH in $(seq -1.5 0.1 1.5)  
do
qsub -N "$STRUC" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,THRESH=$THRESH,NPARC=$NPARC -hold_jid $TP SCSTPtrans.sh
qsub -N "$STRUCNULL" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,THRESH=$THRESH,NPARC=$NPARC -hold_jid $NULLSC nullSCSTPtrans.sh
qsub -N "$BCTSTRUCNULL" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,THRESH=$THRESH,NPARC=$NPARC -hold_jid $NULLSC BCTnullSCSTPtrans.sh

qsub -N "plotstruc" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K -hold_jid "$STRUC" plotStrucTP.sh
qsub -N "plotsubjstrucv2" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,THRESH=$THRESH -hold_jid "$STRUC" plotsubjStrucTPv2.sh
qsub -N "plotnullstruc" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K -hold_jid "$STRUCNULL" plotNULLStrucTP.sh
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K -hold_jid "$BCTSTRUCNULL" plotBCTNULLStrucTP.sh
done

qsub -l h_vmem=8G,s_vmem=7.5G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,K=$K -hold_jid "$STRUC" SCTPthresh.sh

###############
### Control ###
###############

# group average

# subject level

#################
### Cognition ###
#################

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN -hold_jid $TP dprimeTPDurnew.sh
