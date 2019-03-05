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
TP={}
ASSIGN={}
SYMM='symm'

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN systems_plot.sh

#HCP=$(qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,METHOD=$METHOD,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $ASSIGN hcp_cluster.sh)
#HCP="${HCP//[!0-9]/}"
HCP={}
#qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $HCP hcp_tpmat.sh
qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $HCP hcp_dwelltime.sh

BEGINCOMMENT
NULL=$(qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,N=$NPERMS,BD=$BASEDIR,MP=$MATPATH -hold_jid $TP nullTransProbs.sh)
NULL="${NULL//[!0-9]/}"

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid $NULL restvsnbackNP.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP plotDur.sh

qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid $TP nbackblockDur.sh
NBACK="${NBACK//[!0-9]/}"


###################
### Development ###
###################
ENDCOMMENT

#qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid $TP ageTPDur.sh
#qsub -l h_vmem=15.5G,s_vmem=15G -q all.q,basic.q,all.short.q,himem.q -v D=$ROOT,K=$K,BD=$BASEDIR,RP=$RPATH -hold_jid "$SYMM" ageTPprop.sh

##########################
### Structure-function ###
##########################

BEGINCOMMENT
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