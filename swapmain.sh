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
LAB='final'

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

# PROC=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q -v ZDIM=$ZDIM,NPARC=$NPARC,SCAN=$SCAN,LAB=$LAB,BD=$BASEDIR,MP=$MATPATH ProcessDataBaumSample.sh)
# PROC="${PROC//[!0-9]/}"
# NULLSC=$(qsub -l h_vmem=30.5G,s_vmem=30G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,NPARC=$NPARC,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC nullSC.sh)
# NULLSC="${NULLSC//[!0-9]/}"

##########################
### k-means clustering ###
##########################

NREPS=1
NSPLITS=100
REPK='kmeans'$ROOT
SPLITHALVES='splithalves'$ROOT

# for K in {2..11}
# do
# 	for S in $(seq $NSPLITS)
# 	do
# 	qsub -N "$REPK" -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,S=$S,D=$ROOT,N=$NREPS,Z=$ZDIM,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC kmeans.sh
# 	done
# done

##########################
### Process clustering ###
##########################

PROC={}
NULLSC={}
REPK={}

K=5

ASSIGN0=$(qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,qlogin.q,reboot.q -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,BD=$BASEDIR,MP=$MATPATH -hold_jid "$REPK" getassignmentspy.sh)
ASSIGN0="${ASSIGN0//[!0-9]/}"
ASSIGN=$(qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,qlogin.q,reboot.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "$ASSIGN0" reorderClusters.sh)
ASSIGN="${ASSIGN//[!0-9]/}"

###########################################################
### Supplemental control analyses related to clustering ###
###########################################################

qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN kmeansmotion.sh
qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN kmeanstaskregress.sh
qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN distancetostate.sh
qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN inter_RSN_fc.sh

#qsub -N "nullts" -l h_vmem=24.5G,s_vmem=24G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC makenullts.sh
#qsub -N "nullcluster" -l h_vmem=24.5G,s_vmem=24G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "nullts" clusternullts.sh
#qsub -l h_vmem=24.5G,s_vmem=24G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "nullcluster" plot_silhouette_null.sh

#########################
### Assess clustering ###
#########################

for K in {2..11}
do
	qsub -N "zrand" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN calczrand.sh
done

qsub -N "plotzrand" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,BD=$BASEDIR,MP=$MATPATH -hold_jid "zrand" plotzrand.sh

qsub -l h_vmem=10.5G,s_vmem=10G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN stateRepresentation.sh

K=5
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN systems_plot.sh

######################
### State Dynamics ###
######################

NPERMS=5000
TP='transprobs'
SYMM='symm'

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN getcentroids.sh
TP=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN getDynamics.sh)
TP="${TP//[!0-9]/}"

NULL=$(qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,N=$NPERMS,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP nullTransProbs.sh)
NULL="${NULL//[!0-9]/}"
qsub -N "$SYMM" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP symmRvNv2.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP persistNull.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $NULL restvsnbackNP.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP plotDur.sh

NB=$(qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP nbackblockDur.sh)
NB="${NB//[!0-9]/}"
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP TPDistance.sh

###################
### Development ###
###################

qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $NB ageTPDur.sh
qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$SYMM" ageTPprop.sh

##########################
### Structure-function ###
##########################

STRUC='struc'$ROOT
STRUCNULL='strucnull'$ROOT
BCTSTRUCNULL='BCTstrucnull'$ROOT

for THRESH in $(seq -1.5 0.1 1.5)  
do
qsub -N "$STRUC" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP SCSTPtrans.sh
qsub -N "$STRUCNULL" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,MP=$MATPATH -hold_jid $NULLSC,$TP nullSCSTPtrans.sh
qsub -N "$BCTSTRUCNULL" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,MP=$MATPATH -hold_jid $NULLSC,$TP BCTnullSCSTPtrans.sh

qsub -N "plotstruc" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$STRUC" plotStrucTP.sh
qsub -N "plotsubjstrucv2" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,RP=$RPATH -hold_jid "$STRUC" plotsubjStrucTPv2.sh
qsub -N "plotnullstruc" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$STRUCNULL" plotNULLStrucTP.sh
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$BCTSTRUCNULL" plotBCTNULLStrucTP.sh
done

qsub -l h_vmem=8G,s_vmem=7.5G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$BCTSTRUCNULL" SCTPthresh.sh

###############
### Control ###
###############

# group average

DMAT=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC getDmat.sh)
DMAT="${DMAT//[!0-9]/}"

GROUPSC=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $DMAT makeGroupSC.sh)
GROUPSC="${GROUPSC//[!0-9]/}"

NSPLITS_SPHERE=25
#BEGINCOMMENT
for Ki in $(seq $K); do  # loop through clusters and compute spatial correlation-preserving null activity patterns
	for S in $(seq $NSPLITS_SPHERE); do 	# split up computation by parallelizing across jobs
		qsub -N "SPHERE_NULL$S$Ki" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,Ki=$Ki,S=$S,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN make_nullsphere_states.sh
	done
done
#ENDCOMMENT
NSPLITS_DLWNULL=100
NPERMS_DLWNULL=1000
BEGINCOMMENT
for S in $(seq $NSPLITS_DLWNULL); do 	# get distribution of null models for group average SC
	qsub -N "DLW_NULL$S" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_DLWNULL,NPERMS=$NPERMS_DLWNULL,S=$S -hold_jid $DMAT,$GROUPSC precomputeDLWNull.sh
done
ENDCOMMENT

NSPLITS_BOOTSC=100
NPERMS_BOOTSC=1000
BEGINCOMMENT
for S in $(seq $NSPLITS_BOOTSC); do 	# bootstrap group average SC
	qsub -N "BOOTSC$S" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_BOOTSC,NPERMS=$NPERMS_BOOTSC,S=$S -hold_jid $DMAT precomputeBootstrapSC.sh
done
ENDCOMMENT

K=5
# compute persistence energy for distribution of null activity patterns in single null networks
GRP_ENERGY=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,NSPLITS=$NSPLITS_SPHERE,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "SPHERE_NULL$NSPLITS_SPHERE$K",$GROUPSC persistEnergyGroup.sh)
GRP_ENERGY="${GRP_ENERGY//[!0-9]/}"
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid $GRP_ENERGY plotEnergySphereGroup.sh

# compute persistence energy for actual activity patterns in distribution of null networks vs. single group representative network
DLWNULL_ENERGY=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,NPERMS=$NPERMS_DLWNULL,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "DLW_NULL$NSPLITS_DLWNULL" persistEnergyDLWNulls.sh)
DLWNULL_ENERGY="${DLWNULL_ENERGY//[!0-9]/}"
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid $DLWNULL_ENERGY,$GRP_ENERGY plotEnergyDynamicsGroup.sh

# compute correlation between persistence energy and bootstrapped group representative networks
BOOTSC_ENERGY=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,NPERMS=$NPERMS_BOOTSC,K=$K,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "BOOTSC$NSPLITS_BOOTSC" persistEnergyBootstrapGroup.sh)
BOOTSC_ENERGY="${BOOTSC_ENERGY//[!0-9]/}"

#################
### Cognition ###
#################

qsub -l h_vmem=6.5G,s_vmem=6G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $NB dprimeTPDur.sh

###########
### HCP ###
###########

# replicate results in HCP dataset with equal amounts of resting state and n-back data
# must have concatenated HCP data in data folder (see load_hcpdata.m)
BEGINCOMMENT
HCP=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,METHOD=$METHOD,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN hcp_cluster.sh)
HCP="${HCP//[!0-9]/}"
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $HCP hcp_tpmat.sh
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP hcp_dwelltime.sh
ENDCOMMENT