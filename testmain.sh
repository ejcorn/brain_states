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
LAB='troubleshoot'

ROOT='Scan'$SCAN'Laus'$NPARC'Z'$ZDIM$LAB
BASEDIR='/data/tesla-data/ecornblath/brain_states/'
MATPATH=/share/apps/matlab/R2017a/bin/matlab
RPATH=/share/apps/R/R-3.2.5/bin/Rscript
MASTERDIR=$BASEDIR'results/'$ROOT
if [ ! -d "$MASTERDIR" ]; then
  mkdir -p $MASTERDIR		# recursively create general results folder and output folder
fi

cd $BASEDIR'jobs'		# change to directory containing all shell scripts

BEGINCOMMENT
NSPLITS=10
NREPS=1
REPK={}
K=5

ASSIGN0=$(qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,BD=$BASEDIR,MP=$MATPATH -hold_jid "$REPK" getassignmentspy.sh)
ASSIGN0="${ASSIGN0//[!0-9]/}"
ASSIGN=$(qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "$ASSIGN0" reorderClusters.sh)
ASSIGN="${ASSIGN//[!0-9]/}"
ENDCOMMENT

K=5
NPERMS=5000
ASSIGN={}
PROC={}
SYMM='symm'
TP={}
#qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP plotDur.sh

#qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP nbackblockDur.sh
#qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN distancetostate.sh

#qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN inter_RSN_fc.sh

#qsub -N "nullts" -l h_vmem=24.5G,s_vmem=24G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC makenullts.sh
#qsub -N "nullcluster" -l h_vmem=24.5G,s_vmem=24G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "nullts" clusternullts.sh
#qsub -l h_vmem=24.5G,s_vmem=24G -q all.q,basic.q,all.short.q,himem.q -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "nullcluster" plot_silhouette_null.sh

#TP=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN getDynamics.sh)
#TP="${TP//[!0-9]/}"

###############
### Control ###
###############

K=5
BEGINCOMMENT
DMAT=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $PROC getDmat.sh)
DMAT="${DMAT//[!0-9]/}"
DMAT={}

GROUPSC=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid $DMAT makeGroupSC.sh)
GROUPSC="${GROUPSC//[!0-9]/}"

NSPLITS_SPHERE=25
for Ki in $(seq $K); do
	for S in $(seq $NSPLITS_SPHERE); do
		SPHERE_NULL=$(qsub -l h_vmem=16.5G,s_vmem=16G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,Ki=$Ki,S=$S,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN make_nullsphere_states.sh)
		SPHERE_NULL="${SPHERE_NULL//[!0-9]/}"
	done
done


NSPLITS_GRAM=50
for S in $(seq $NSPLITS_GRAM); do
	qsub -N "Gramian$S" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_GRAM,S=$S -hold_jid $PROC precomputeGramian.sh
done


ENDCOMMENT
NSPLITS_DLWNULL=100
NPERMS_DLWNULL=1000
BEGINCOMMENT
for S in $(seq $NSPLITS_DLWNULL); do
	qsub -N "DLW_NULL$S" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_DLWNULL,NPERMS=$NPERMS_DLWNULL,S=$S -hold_jid $DMAT precomputeDLWNull.sh
done
ENDCOMMENT
BEGINCOMMENT
NSPLITS_BOOTSC=100
NPERMS_BOOTSC=1000
for S in $(seq $NSPLITS_BOOTSC); do
	qsub -N "BOOTSC$S" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_BOOTSC,NPERMS=$NPERMS_BOOTSC,S=$S -hold_jid $DMAT precomputeBootstrapSC.sh
done

ENDCOMMENT
NSPLITS_SPHERE=25
SPHERE_NULL={}
GROUPSC={}

GRP_ENERGY=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,NSPLITS=$NSPLITS_SPHERE,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $SPHERE_NULL,$GROUPSC persistEnergyGroup.sh)
GRP_ENERGY="${GRP_ENERGY//[!0-9]/}"

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid $GRP_ENERGY plotEnergySphereGroup.sh

NSPLITS_GRAM=50

#SUBJ_ENERGY=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "Gramian$NSPLITS_GRAM" persistEnergySubjects.sh)
#SUBJ_ENERGY="${SUBJ_ENERGY//[!0-9]/}"

DLW_NULL={}
DLWNULL_ENERGY=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,NPERMS=$NPERMS_DLWNULL,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "$DLW_NULL$NSPLITS_DLWNULL" persistEnergyDLWNulls.sh)
DLWNULL_ENERGY="${DLWNULL_ENERGY//[!0-9]/}"
#DLWNULL_ENERGY={}

NSPLITS_BOOTSC=100

NPERMS_BOOTSC=1000
#BOOTSC={}
#BOOTSC_ENERGY=$(qsub -l h_vmem=20.5G,s_vmem=20G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,NPERMS=$NPERMS_BOOTSC,K=$K,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "BOOTSC$NSPLITS_BOOTSC" persistEnergyBootstrapGroup.sh)
#BOOTSC_ENERGY="${BOOTSC_ENERGY//[!0-9]/}"

DLWNULL_ENERGY={}

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid $DLWNULL_ENERGY,$GRP_ENERGY plotEnergyDynamicsGroup.sh

#################
### Cognition ###
#################
BEGINCOMMENT
#NB=$(qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP nbackblockDur.sh)
#NB="${NB//[!0-9]/}"
NB={}
qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP ageTPDur.sh
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $NB dprimeTPDur.sh

##########################
### Structure-function ###
##########################


STRUC='struc'$ROOT
STRUCNULL='strucnull'$ROOT
BCTSTRUCNULL='BCTstrucnull'$ROOT
NULLSC={}

for THRESH in $(seq -1.5 0.1 1.5)  
do
	qsub -N "$STRUC" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP SCSTPtrans.sh
	qsub -N "$BCTSTRUCNULL" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,MP=$MATPATH -hold_jid $NULLSC BCTnullSCSTPtrans.sh
	qsub -N "plotstruc" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$STRUC" plotStrucTP.sh
	qsub -N "plotsubjstrucv2" -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,THRESH=$THRESH,BD=$BASEDIR,RP=$RPATH -hold_jid "$STRUC" plotsubjStrucTPv2.sh
	qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,THRESH=$THRESH,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$BCTSTRUCNULL" plotBCTNULLStrucTP.sh
done

qsub -l h_vmem=8G,s_vmem=7.5G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,SCAN=$SCAN,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$STRUC" SCTPthresh.sh
ENDCOMMENT