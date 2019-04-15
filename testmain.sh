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

# PROC=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q -v ZDIM=$ZDIM,NPARC=$NPARC,SCAN=$SCAN,LAB=$LAB,BD=$BASEDIR,MP=$MATPATH ProcessDataBaumSample.sh)
# PROC="${PROC//[!0-9]/}"
# NULLSC=$(qsub -l h_vmem=30.5G,s_vmem=30G -q $QUEUE -v D=$ROOT,SCAN=$SCAN,NPARC=$NPARC,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC nullSC.sh)
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
# 	qsub -N "$REPK" -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,S=$S,D=$ROOT,N=$NREPS,Z=$ZDIM,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC kmeans.sh
# 	done
# done

##########################
### Process clustering ###
##########################


NULLSC={}
REPK={}
PROC={}

K=5

GROUPSC={}
DMAT={}

NSPLITS_SPHERE=50
BEGINCOMMENT
for Ki in $(seq $K); do  # loop through clusters and compute spatial correlation-preserving null activity patterns
	for S in $(seq $NSPLITS_SPHERE); do 	# split up computation by parallelizing across jobs
		qsub -N "SPHERE_NULL$S$Ki" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,Ki=$Ki,S=$S,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN make_nullsphere_states.sh
	done
done
ENDCOMMENT
NSPLITS_DLWNULL=100
NPERMS_DLWNULL=1000
BEGINCOMMENT
for S in $(seq $NSPLITS_DLWNULL); do 	# get distribution of null models for group average SC
	qsub -N "DLW_NULL$S" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_DLWNULL,NPERMS=$NPERMS_DLWNULL,S=$S -hold_jid $DMAT,$GROUPSC precomputeDLWNull.sh
done
ENDCOMMENT

NSPLITS_BOOTSC=100
NPERMS_BOOTSC=1000
BEGINCOMMENT
for S in $(seq $NSPLITS_BOOTSC); do 	# bootstrap group average SC
	qsub -N "BOOTSC$S" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_BOOTSC,NPERMS=$NPERMS_BOOTSC,S=$S -hold_jid $DMAT precomputeBootstrapSC.sh
done
ENDCOMMENT

NSPLITS_SWEEP=25
for S in $(seq $NSPLITS_SWEEP); do 	# sweep through control horizon and normalization for each bootstrapped SC matrix
	qsub -N "SWEEP$S" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,K=$K,MP=$MATPATH,NSPLITS=$NSPLITS_SWEEP,NPERMS=$NPERMS_BOOTSC,S=$S -hold_jid $DMAT,"BOOTSC$NSPLITS_BOOTSC" scalinghorizonsweep.sh
done
qsub -N "SWEEPCORR" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NSPLITS=$NSPLITS_SWEEP,NPERMS=$NPERMS_BOOTSC,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "SWEEP$NSPLITS_SWEEP" sweepcorr.sh
qsub -l h_vmem=3.5G,s_vmem=3G -q $QUEUE -v D=$ROOT,NPERMS=$NPERMS_BOOTSC,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "SWEEPCORR" plotEnergyDynamicsSweep.sh

K=5

# loop through a few control horizons and normalizations to make sure results are consistent
T_RNG=(1 5 10)
for c in $(seq 0 5 10); do
	for T in "${T_RNG[@]}"; do
		# compute persistence energy for distribution of null activity patterns in single null networks
		qsub -N "PESPHERE$c$T" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NSPLITS=$NSPLITS_SPHERE,K=$K,c=$c,T=$T,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "SPHERE_NULL$NSPLITS_SPHERE$K",$GROUPSC persistEnergyGroup.sh
	done
done

c=0
for T in "${T_RNG[@]}"; do
	# compute persistence energy for actual activity patterns in distribution of null networks vs. single group representative network
	# can only use c = 0 b/c that's what i precomputed and you need to normalize before making null
	qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NPERMS=$NPERMS_DLWNULL,K=$K,c=$c,T=$T,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "DLW_NULL$NSPLITS_DLWNULL","PESPHERE$c$T" persistEnergyDLWNulls.sh
done
