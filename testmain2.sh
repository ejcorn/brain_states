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
QUEUE=all.q,basic.q,all.short.q,himem.q,qlogin.q,qlogin.long.q,qlogin.himem.q
ROOT='Scan'$SCAN'Laus'$NPARC'Z'$ZDIM$LAB
BASEDIR='/data/tesla-data/ecornblath/brain_states/'
MATPATH=/share/apps/matlab/R2017a/bin/matlab
RPATH=/share/apps/R/R-3.2.5/bin/Rscript
MASTERDIR=$BASEDIR'results/'$ROOT
if [ ! -d "$MASTERDIR" ]; then
  mkdir -p $MASTERDIR		# recursively create general results folder and output folder
fi

cd $BASEDIR'jobs'		# change to directory containing all shell scripts

NREPS=1
NSPLITS=10
REPK={}
K=5

################
### Analysis ###
################

NPERMS=5000
TP='transprobs'
SYMM='symm'
K=5

##########
### TP ###
##########

ASSIGN={}
NULL={}
TP={}
NB={}

#qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP persistNull.sh
#qsub -l h_vmem=6.5G,s_vmem=6G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $NB dprimeTPDur.sh
#qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $NB ageTPDur.sh
#qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$SYMM" ageTPprop.sh



###############
### Control ###
###############


# group average
BEGINCOMMENT
DMAT=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC getDmat.sh)
DMAT="${DMAT//[!0-9]/}"

GROUPSC=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $DMAT makeGroupSC.sh)
GROUPSC="${GROUPSC//[!0-9]/}"
ENDCOMMENT
GROUPSC={}
NSPLITS_SPHERE=25
BEGINCOMMENT
for Ki in $(seq $K); do  # loop through clusters and compute spatial correlation-preserving null activity patterns
	for S in $(seq $NSPLITS_SPHERE); do 	# split up computation by parallelizing across jobs
		qsub -N "SPHERE_NULL$S$Ki" -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,Ki=$Ki,S=$S,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN make_nullsphere_states.sh
	done
done
ENDCOMMENT
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
GRP_ENERGY=$(qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NSPLITS=$NSPLITS_SPHERE,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "SPHERE_NULL$NSPLITS_SPHERE$K",$GROUPSC,5065512 persistEnergyGroup.sh)
GRP_ENERGY="${GRP_ENERGY//[!0-9]/}"
qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid $GRP_ENERGY plotEnergySphereGroup.sh
BEGINCOMMENT
# compute persistence energy for actual activity patterns in distribution of null networks vs. single group representative network
DLWNULL_ENERGY=$(qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NPERMS=$NPERMS_DLWNULL,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "DLW_NULL$NSPLITS_DLWNULL",5065512 persistEnergyDLWNulls.sh)
DLWNULL_ENERGY="${DLWNULL_ENERGY//[!0-9]/}"
qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid $DLWNULL_ENERGY,$GRP_ENERGY plotEnergyDynamicsGroup.sh

# compute correlation between persistence energy and bootstrapped group representative networks
BOOTSC_ENERGY=$(qsub -l h_vmem=20.5G,s_vmem=20G -q $QUEUE -v D=$ROOT,NPERMS=$NPERMS_BOOTSC,K=$K,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "BOOTSC$NSPLITS_BOOTSC",5065512 persistEnergyBootstrapGroup.sh)
BOOTSC_ENERGY="${BOOTSC_ENERGY//[!0-9]/}"
ENDCOMMENT