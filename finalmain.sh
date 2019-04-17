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
METHOD='cosine'
NPARC=250
SCAN='C'
LAB='cosinefinal'

ROOT='Scan'$SCAN'Laus'$NPARC'Z'$ZDIM$LAB
BASEDIR='/data/tesla-data/ecornblath/brain_states/'
MATPATH=/share/apps/matlab/R2017a/bin/matlab
RPATH=/share/apps/R/R-3.2.5/bin/Rscript
QUEUE=all.q,basic.q,all.short.q,himem.q,qlogin.q,qlogin.long.q,qlogin.himem.q
MASTERDIR=$BASEDIR'results/'$ROOT
if [ ! -d "$MASTERDIR" ]; then
  mkdir -p $MASTERDIR		# recursively create general results folder and output folder
fi

cd $BASEDIR'jobs'		# change to directory containing all shell scripts

####################
### Process data ###
####################

PROC=$(qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v ZDIM=$ZDIM,NPARC=$NPARC,SCAN=$SCAN,LAB=$LAB,BD=$BASEDIR,MP=$MATPATH ProcessDataBaumSample.sh)
PROC="${PROC//[!0-9]/}"
NULLSC=$(qsub -l h_vmem=30.5G,s_vmem=30G -q $QUEUE -v D=$ROOT,SCAN=$SCAN,NPARC=$NPARC,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC nullSC.sh)
NULLSC="${NULLSC//[!0-9]/}"

##########################
### k-means clustering ###
##########################

NREPS=1
NSPLITS=20
REPK='kmeans'$ROOT
SPLITHALVES='splithalves'$ROOT

for K in {2..11}
do
	for S in $(seq $NSPLITS)
	do
	qsub -N "$REPK" -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,S=$S,D=$ROOT,N=$NREPS,Z=$ZDIM,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC kmeans.sh
	done
done

##########################
### Process clustering ###
##########################

K=5

ASSIGN0=$(qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,BD=$BASEDIR,MP=$MATPATH -hold_jid "$REPK" getassignmentspy.sh)
ASSIGN0="${ASSIGN0//[!0-9]/}"
ASSIGN=$(qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "$ASSIGN0" reorderClusters.sh)
ASSIGN="${ASSIGN//[!0-9]/}"

###########################################################
### Supplemental control analyses related to clustering ###
###########################################################

qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN kmeansmotion.sh
qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN kmeanstaskregress.sh
qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN distancetostate.sh
qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN inter_RSN_fc.sh

for K in {2..11}; do
	qsub -N "DSTASK$K" -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN dstaskcluster.sh
done
K=5
qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "DSTASK*" dstaskassess.sh

NSPLITS_SH=50
for S in $(seq 1 $NSPLITS_SH); do
	qsub -N "SH$S" -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,S=$S,D=$ROOT,Z=$ZDIM,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC SHkmeans.sh
done
qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,NSPLITS=$NSPLITS_SH,BD=$BASEDIR,MP=$MATPATH -hold_jid "SH*" analyzeSHkmeans.sh

#qsub -N "nullts" -l h_vmem=24.5G,s_vmem=24G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC makenullts.sh
#qsub -N "nullcluster" -l h_vmem=24.5G,s_vmem=24G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "nullts" clusternullts.sh
#qsub -l h_vmem=24.5G,s_vmem=24G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "nullcluster" plot_silhouette_null.sh

#########################
### Assess clustering ###
#########################

for K in {2..11}
do
	qsub -N "zrand" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN calczrand.sh
done

qsub -N "plotzrand" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,BD=$BASEDIR,MP=$MATPATH -hold_jid "zrand" plotzrand.sh

qsub -l h_vmem=10.5G,s_vmem=10G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN stateRepresentation.sh

K=5
qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN systems_plot.sh

######################
### State Dynamics ###
######################

NPERMS=5000
TP='transprobs'
SYMM='symm'

qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN getcentroids.sh
TP=$(qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN getDynamics.sh)
TP="${TP//[!0-9]/}"

NULL=$(qsubp4 -l h_vmem=15.5G,s_vmem=15G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,N=$NPERMS,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP nullTransProbs.sh)
NULL="${NULL//[!0-9]/}"
qsub -N "$SYMM" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP symmRvNv2.sh
qsub -N "$SYMM" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $ASSIGN RvNMI.sh

qsubp4 -l h_vmem=20.5G,s_vmem=20G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP persistNull.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $NULL restvsnbackNP.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP plotDur.sh

NB=$(qsub -l h_vmem=15.5G,s_vmem=15G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP nbackblockDur.sh)
NB="${NB//[!0-9]/}"
qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP TPDistance.sh

###################
### Development ###
###################

qsub -l h_vmem=15.5G,s_vmem=15G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $NB ageTPDur.sh
qsub -l h_vmem=15.5G,s_vmem=15G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$SYMM" ageTPprop.sh

##########################
### Structure-function ###
##########################

STRUC='struc'$ROOT
STRUCNULL='strucnull'$ROOT
BCTSTRUCNULL='BCTstrucnull'$ROOT

for THRESH in $(seq -1.5 0.1 1.5)  
do
qsub -N "$STRUC" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP SCSTPtrans.sh
qsub -N "$STRUCNULL" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,MP=$MATPATH -hold_jid $NULLSC,$TP nullSCSTPtrans.sh
qsub -N "$BCTSTRUCNULL" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,MP=$MATPATH -hold_jid $NULLSC,$TP BCTnullSCSTPtrans.sh

qsub -N "plotstruc" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$STRUC" plotStrucTP.sh
qsub -N "plotsubjstrucv2" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,RP=$RPATH -hold_jid "$STRUC" plotsubjStrucTPv2.sh
qsub -N "plotnullstruc" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$STRUCNULL" plotNULLStrucTP.sh
qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$BCTSTRUCNULL" plotBCTNULLStrucTP.sh
done

qsub -l h_vmem=8G,s_vmem=7.5G -q $QUEUE -v D=$ROOT,SCAN=$SCAN,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$BCTSTRUCNULL" SCTPthresh.sh

###############
### Control ###
###############

# group average

DMAT=$(qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN getDmat.sh)
DMAT="${DMAT//[!0-9]/}"

GROUPSC=$(qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $DMAT makeGroupSC.sh)
GROUPSC="${GROUPSC//[!0-9]/}"

NSPLITS_SPHERE=50
#BEGINCOMMENT
for Ki in $(seq $K); do  # loop through clusters and compute spatial correlation-preserving null activity patterns
	for S in $(seq $NSPLITS_SPHERE); do 	# split up computation by parallelizing across jobs
		qsub -N "SPHERE_NULL$S$Ki" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,Ki=$Ki,S=$S,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN make_nullsphere_states.sh
	done
done
#ENDCOMMENT
NSPLITS_DLWNULL=100
NPERMS_DLWNULL=1000
#BEGINCOMMENT
for S in $(seq $NSPLITS_DLWNULL); do 	# get distribution of null models for group average SC
	qsub -N "DLW_NULL$S" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_DLWNULL,NPERMS=$NPERMS_DLWNULL,S=$S -hold_jid $DMAT,$GROUPSC precomputeDLWNull.sh
done
#ENDCOMMENT

NSPLITS_BOOTSC=100
NPERMS_BOOTSC=1000
#BEGINCOMMENT
for S in $(seq $NSPLITS_BOOTSC); do 	# bootstrap group average SC
	qsub -N "BOOTSC$S" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_BOOTSC,NPERMS=$NPERMS_BOOTSC,S=$S -hold_jid $DMAT precomputeBootstrapSC.sh
done
#ENDCOMMENT

for S in $(seq $NSPLITS_SWEEP); do 	# sweep through control horizon and normalization for each bootstrapped SC matrix
	qsub -N "SWEEP$S" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,K=$K,MP=$MATPATH,NSPLITS=$NSPLITS_SWEEP,NPERMS=$NPERMS_BOOTSC,S=$S -hold_jid $DMAT,"BOOTSC*" scalinghorizonsweep.sh
done

qsub -N "SWEEPCORR" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NSPLITS=$NSPLITS_SWEEP,NPERMS=$NPERMS_BOOTSC,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "SWEEP*" sweepcorr.sh
qsub -l h_vmem=3.5G,s_vmem=3G -q $QUEUE -v D=$ROOT,NPERMS=$NPERMS_BOOTSC,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "SWEEPCORR" plotEnergyDynamicsSweep.sh

# loop through a few control horizons and normalizations to make sure results are consistent
T_RNG=(1 5 10)
C_RNG=(0 1 5 10)
for c in "${C_RNG[@]}"; do
	for T in "${T_RNG[@]}"; do
		# compute persistence energy for distribution of null activity patterns in single null networks
		qsub -N "PESPHERE$c$T" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NSPLITS=$NSPLITS_SPHERE,K=$K,c=$c,T=$T,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "SPHERE_NULL*",$GROUPSC persistEnergyGroup.sh
	done
done

c=0
for T in "${T_RNG[@]}"; do
	# compute persistence energy for actual activity patterns in distribution of null networks vs. single group representative network
	# can only use c = 0 b/c that's what i precomputed and you need to normalize before making null
	qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NPERMS=$NPERMS_DLWNULL,K=$K,c=$c,T=$T,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "DLW_NULL*","PESPHERE*" persistEnergyDLWNulls.sh
done

#################
### Cognition ###
#################

qsub -l h_vmem=6.5G,s_vmem=6G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $NB dprimeTPDur.sh

###########
### HCP ###
###########

# replicate results in HCP dataset with equal amounts of resting state and n-back data
# must have concatenated HCP data in data folder (see load_hcpdata.m)
BEGINCOMMENT
HCP=$(qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,METHOD=$METHOD,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN hcp_cluster.sh)
HCP="${HCP//[!0-9]/}"
qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $HCP hcp_tpmat.sh
qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP hcp_dwelltime.sh
ENDCOMMENT