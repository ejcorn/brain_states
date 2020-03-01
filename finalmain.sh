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
NPARC=250 # to reproduce Fig S5, change this to 125 and run everything below
SCAN='C'
LAB='prepub_test' # this can be anything, just an identifier for output folder

ROOT='Scan'$SCAN'Laus'$NPARC'Z'$ZDIM$LAB
BASEDIR='/data/tesla-data/ecornblath/brain_states_prepubtest/'
MATPATH=/share/apps/matlab/R2017a/bin/matlab
MATPATH2014=matlab # use earlier version of matlab to use less memory and get slots
RPATH=/share/apps/R/R-3.2.5/bin/Rscript
QUEUE=all.q,basic.q,all.short.q,himem.q,qlogin.q,qlogin.long.q,qlogin.himem.q
LONGJOB=h_rt=99:00:00,s_rt=99:00:00 # set long time limit for CPU dump if job actually takes a while
MASTERDIR=$BASEDIR'results/'$ROOT
if [ ! -d "$MASTERDIR" ]; then
  mkdir -p $MASTERDIR		# recursively create general results folder and output folder
fi

cd $BASEDIR'jobs'		# change to directory containing all shell scripts

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~ MAIN ~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# carry out analyses in main text of paper

####################
### Process data ###
####################

# load preprocessed PNC BOLD and DTI data

qsub -N "PROC" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v ZDIM=$ZDIM,NPARC=$NPARC,SCAN=$SCAN,LAB=$LAB,BD=$BASEDIR,MP=$MATPATH ProcessDataBaumSample.sh

##########################
### k-means clustering ###
##########################
# perform many parallel iterations of k-means clustering for k = 2 to 11

NREPS=1
NSPLITS=20

for K in {2..11}
do
	for S in $(seq $NSPLITS)
	do
	qsub -N "REPK${K}${S}" -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,S=$S,D=$ROOT,N=$NREPS,Z=$ZDIM,BD=$BASEDIR,MP=$MATPATH -hold_jid "PROC" kmeans.sh
	done
done

##########################
### Process clustering ###
##########################

# process clustering output: concatenate repetitions of clustering, reorder cluster indices from arbitrary to predefined a set order for consistent visualization at k = 5

K=5 # to reproduce Fig S12, just change this to 6 and run everything below

qsub -N "ASSIGNCLUSTERS_INIT" -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,D=$ROOT,N=$NSPLITS,BD=$BASEDIR,MP=$MATPATH -hold_jid "REPK*" getassignmentspy.sh
qsub -N "ASSIGNCLUSTERS_FINAL" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "ASSIGNCLUSTERS_INIT" reorderClusters.sh

#########################
### Assess clustering ###
#########################

# evaluate number of clusters using elbow criterion and change in variance explained per Gutierrez-Barragan et al. 2019. --> Fig S2a-b
qsub -N "elbow" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v METHOD=$METHOD,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "ASSIGNCLUSTERS_FINAL" elbow.sh

# compare representation of clusters between rest and n-back, Fig S2
qsub -l h_vmem=10.5G,s_vmem=10G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,METHOD=$METHOD -hold_jid "ASSIGNCLUSTERS_FINAL" stateRepresentation.sh

# compare alignment of cluster centroids with Yeo resting state networks --> Fig 2b
qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "ASSIGNCLUSTERS_FINAL" systems_plot.sh

######################
### State Dynamics ###
######################

# name and save cluster centroids for visualization using localplotcentroids.sh --> Fig 2a, Fig. S4a, Fig. S4b
# these cluster centroids will be used for multiple analyses throughout, expecially control theory in Fig. 5
qsub -N "CENTROIDS$K" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH -hold_jid "ASSIGNCLUSTERS_FINAL" getcentroids.sh

# compute metrics of brain state dynamics --> Fig 3-4
qsub -N "GETDYNAMICS$K" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH -hold_jid "ASSIGNCLUSTERS_FINAL" getDynamics.sh 

qsub -l h_vmem=2.5G,s_vmem=2G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid "GETDYNAMICS$K" plotDur.sh # Fig 3a-c

qsub -N "NBACK$K" -l h_vmem=15.5G,s_vmem=15G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "GETDYNAMICS$K" nbackblockDur.sh # data for Fig 3d

qsub -N "PLOTNBACK$K" -l h_vmem=15.5G,s_vmem=15G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "NBACK$K" plotnbackblockDur.sh # plots for Fig 3d

qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "NBACK$K" restvs2backNP.sh # Fig 4a-c

#################
### Cognition ###
#################

# relationship between fractional occupancy and WM performance --> Fig 3e
# relationship between transition probabilities and WM performance --> Fig 4d

qsub -l h_vmem=6.5G,s_vmem=6G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH,MP=$MATPATH -hold_jid "NBACK$K" dprimeTPDur.sh

###################
### Development ###
###################

# relationship between age and fractional occupancy, dwell time, transition probabilities --> Fig 6a,b,d-h
qsub -l h_vmem=3.5G,s_vmem=3G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid "NBACK$K" ageTPDur.sh

######################
### Control theory ###
######################

# assess whether control properties of white matter networks explain empirically observed brain dynamics

# compute ROIxROI distance matrix then get group-representative SC using distance-dependendent consistency thresholding
# this group-representative SC matrix is used as the A matrix in equation (1) of main text and supplement
qsub -N "DMAT" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "PROC" getDmat.sh
qsub -N "GROUPSC" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "DMAT" makeGroupSC.sh

# get distribution of geometry + topology-preserving null models (SP null, see methods section and results section "Control properties of white matter networks explain brain state transitions") 
# for group average SC (NSPLITS is how many per job, NPERMS is how many total null matrices)

NSPLITS_DLWNULL=100
NPERMS_DLWNULL=1000

for S in $(seq $NSPLITS_DLWNULL); do
	echo $S
	qsub -N "DLW_NULL$S" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,NSPLITS=$NSPLITS_DLWNULL,NPERMS=$NPERMS_DLWNULL,S=$S -hold_jid "DMAT","GROUPSC" precomputeDLWNull.sh
done

# sweep through values of T to identify time scale at which transition energy best explains transition probabilities
C_RNG=(0) # method of stabilization of A matrix. 0 is most commonly used for continuous LTI systems but see Karrer et al. 2019 (https://arxiv.org/abs/1908.03514) for discussion
T_RNG=(5) # setting this to 5 based on the result of Fig S11a... we identified this as T value where structure explains function at rest, so using this T value only for certain analyses to avoid unnecessary computations
SYSTEMS=('VIS' 'SOM' 'DAT' 'VAT' 'LIM' 'FPN' 'DMN')
NPERMS_DLWNULL=1000

for c in "${C_RNG[@]}"; do	
	# compute transition energies for real brain networks and degree-preserving null models across range of T values
	
	qsub -N "TSWEEP$K$c" -l h_vmem=12.5G,s_vmem=12G,$LONGJOB -q $QUEUE -v D=$ROOT,K=$K,c=$c,BD=$BASEDIR,MP=$MATPATH -hold_jid "GROUPSC","CENTROIDS$K" transitionEnergyTSweep.sh
	
	for T in "${T_RNG[@]}"; do	
		# compute transition energies for SP Null matrices at single T-value --> null testing for Fig 5c
		# compare transition energy magnitudes between SP Null, DP Null, and real networks --> Fig 5b heatmap		
		qsub -N "TEDLWNULL$c$T" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NPERMS=$NPERMS_DLWNULL,K=$K,c=$c,T=$T,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "DLW_NULL*","TSWEEP$K$c" transitionEnergyDLWNulls.sh
		
		# compare transition energies for each subject with age --> Fig 6c
		qsub -N "TE_SUBJECTS_FA$K$c$T" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,c=$c,T=$T,MP=$MATPATH,RP=$RPATH -hold_jid "CENTROIDS$K" transitionEnergySubjects_FA.sh
	done
	
	for SYS in {0..6}; do
		# compute weighted control energies for T sweeps in Fig S11s
		qsub -N "WEIGHTED_CONTROL_$SYS" -l h_vmem=12.5G,s_vmem=12G,$LONGJOB -q $QUEUE -v D=$ROOT,K=$K,c=$c,SYS=$SYS,SYSNAME=${SYSTEMS[$SYS]},BD=$BASEDIR,MP=$MATPATH -hold_jid "GROUPSC","CENTROIDS$K" weightedTransitionEnergyTSweep.sh
	done
	
done

c=0 # this should stay as 0
# in the run used for the paper we identified T=5 as the value of T for which the correlation between
# resting state transition probabilities and transition energies was maximal.
# This value is identified in plotEnergyVsTransitionProbability.R using output from transitionEnergyDynamicsGroupv2_TSweep_CFN.m
# upon a second run, this optimal value was 5.5. you can either change T here to whatever value is optimal
# or go into plotEnergyvsTransitionProbability.R and set T.opt.WB = 5 as I have done,
# which will give you similar results even if 5 isn't the exact optimal value

T=5 

# calculate energies in null models for Fig 5d
NSPLITS=50 # parallelize min control energy calculations across jobs
CONTROL='Weighted'

for S in $(seq $NSPLITS); do
	# weighted minimum control energy -- VIS #
	SYSTEM=VIS
	T=0.001 # this T-value identified through Fig S11b, first panel... this is essentially equivalent to computing inter-state distance in state space
	qsub -N "TEDLWNULL_SPLIT$S$SYSTEM$c$T" -l h_vmem=12.5G,s_vmem=12G,$LONGJOB -q $QUEUE -v D=$ROOT,NPERMS=$NPERMS_DLWNULL,NSPLITS=$NSPLITS,S=$S,K=$K,c=$c,T=$T,CONTROL=$CONTROL,SYSTEM=$SYSTEM,BD=$BASEDIR,MP=$MATPATH2014 -hold_jid "DLW_NULL*" transitionEnergySystemWeightedControlDLWNull.sh
	qsub -N "TEBCTNULL_SPLIT$S$SYSTEM$c$T" -l h_vmem=12.5G,s_vmem=12G,$LONGJOB -q $QUEUE -v D=$ROOT,NPERMS=$NPERMS_DLWNULL,NSPLITS=$NSPLITS,S=$S,K=$K,c=$c,T=$T,CONTROL=$CONTROL,SYSTEM=$SYSTEM,BD=$BASEDIR,MP=$MATPATH2014,RP=$RPATH -hold_jid "GROUPSC" transitionEnergySystemWeightedControlBCTNull.sh	
done

# plot uniformly weighted (B is identity matrix) transition energy vs. transition probability --> Fig 5c
qsub -l h_vmem=3.5G,s_vmem=3G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,c=$c,RP=$RPATH -hold_jid "TEDLWNULL$c*","TSWEEP$K$c" plotEvsTP_Uniform.sh 

# plot VIS-weighted transition energy vs. transition probability --> Fig 5d
qsub -l h_vmem=3.5G,s_vmem=3G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,c=$c,NPERMS=$NPERMS_DLWNULL,RP=$RPATH -hold_jid "WEIGHTED_CONTROL_*","TE*NULL_SPLIT*" plotEvsTP_Weighted.sh

# plot results of all T sweeps -> Fig S11a-d
qsub -l h_vmem=3.5G,s_vmem=3G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,c=$c,RP=$RPATH -hold_jid "WEIGHTED_CONTROL_*","TSWEEP$K$c" plotEvsTPSweep_All.sh


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~ SUPPLEMENT ~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# carry out supplemental analyses

###########################################################
### Supplemental control analyses related to clustering ###
###########################################################

# split reliability of clustering solution --> Fig S2c-f
NSPLITS_SH=50
for S in $(seq 1 $NSPLITS_SH); do
	qsub -N "SH$S" -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,S=$S,D=$ROOT,Z=$ZDIM,BD=$BASEDIR,MP=$MATPATH -hold_jid "PROC" SHkmeans.sh
done
qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,NSPLITS=$NSPLITS_SH,BD=$BASEDIR,MP=$MATPATH -hold_jid "SH*" analyzeSHkmeans.sh

# compare brain states to static FC --> Fig. S8
qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "ASSIGNCLUSTERS_FINAL" inter_RSN_fc.sh


# perform clustering on null data --> Fig S3
qsub -N "nullts" -l h_vmem=24.5G,s_vmem=24G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "PROC" makenullts.sh

qsub -N "nullcluster" -l h_vmem=24.5G,s_vmem=24G,$LONGJOB -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "nullts" clusternullts.sh
qsub -N "nullcentroids" -l h_vmem=24.5G,s_vmem=24G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "nullcluster" iprnullcentroids.sh
qsub -l h_vmem=24.5G,s_vmem=24G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "nullcentroids" plot_silhouette_null.sh


# compare similarity between resting state and task centroids --> Fig S4c-f
# dstaskcluster.sh parallelizes clustering
# dstaskassess.sh analyzes and makes figures
for Ki in {2..11}; do
	qsub -N "DSTASK$Ki" -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$Ki,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "ASSIGNCLUSTERS_FINAL" dstaskcluster.sh
done
qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "DSTASK*" dstaskassess.sh


# repeat clustering after removing high motion frames --> Fig S7
qsub -N "motioncluster" -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "ASSIGNCLUSTERS_FINAL" kmeansmotion.sh
qsub -l h_vmem=16.5G,s_vmem=16G -q $QUEUE -v METHOD=$METHOD,K=$K,D=$ROOT,BD=$BASEDIR,MP=$MATPATH -hold_jid "motioncluster" plotkmeansmotion.sh

#################################################################
### Supplemental analyses related to transition probabilities ###
#################################################################

# test whether transition probabilities occur randomly --> Fig S9a-c
NPERMS=5000
qsubp4 -N "NULLTP" -l h_vmem=15.5G,s_vmem=15G,$LONGJOB -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,N=$NPERMS,BD=$BASEDIR,MP=$MATPATH -hold_jid "GETDYNAMICS$K" nullTransProbs.sh

qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "NULLTP" restvsnbackNP.sh # compare rest and n-back transition probabilities along whole task block

qsubp4 -l h_vmem=20.5G,s_vmem=20G,$LONGJOB -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "GETDYNAMICS$K" persistNull.sh

# assess properties of transition matrices and brain state time series --> Fig S9d-f
# specifically: symmetry, dependence on state-space distances, and auto-mutual information

qsub -N "SYMM" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "GETDYNAMICS$K" symmRvNv2.sh
qsub -N "MI" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "ASSIGNCLUSTERS_FINAL" RvNMI.sh
qsub -N "DIST" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,SCAN=$SCAN,BD=$BASEDIR,RP=$RPATH -hold_jid "GETDYNAMICS$K" TPDistance.sh


#######################
### HCP replication ###
#######################

# replicate results in HCP dataset with equal amounts of resting state and n-back data --> Fig S6a-f
# must have concatenated HCP data in data folder (see load_hcpdata.m)
qsub -N "HCP" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,METHOD=$METHOD,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "ASSIGNCLUSTERS_FINAL" hcp_cluster.sh
qsub -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,BD=$BASEDIR,MP=$MATPATH -hold_jid "HCP" hcp_tpmat.sh

###############################################
### Specificity of state persistence energy ###
###############################################

# functional specificity of low energy required to maintain each state --> Fig S10b

NSPLITS_SPHERE=50

for Ki in $(seq $K); do  # loop through clusters and compute spatial correlation-preserving null activity patterns
	for S in $(seq $NSPLITS_SPHERE); do 	# split up computation by parallelizing across jobs
		qsub -N "SPHERE_NULL$S$Ki" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,K=$K,Ki=$Ki,S=$S,BD=$BASEDIR,MP=$MATPATH -hold_jid "CENTROIDS$K","ASSIGNCLUSTERS_FINAL" make_nullsphere_states.sh
	done
done

C_RNG=(0) # method of stabilization of A matrix. 0 is most commonly used for continuous LTI systems but see Karrer et al. 2019 (https://arxiv.org/abs/1908.03514) for discussion
T_RNG=(5) # setting this to 5 based on the result of Fig S11a... we identified this as T value where structure explains function at rest, so using this T value only for certain analyses to avoid unnecessary computations
for c in "${C_RNG[@]}"; do
	for T in "${T_RNG[@]}"; do
		# compute persistence energy for distribution of null activity patterns in single null networks
		# make plot for Figure S10b
		qsub -N "PESPHERE$c$T" -l h_vmem=12.5G,s_vmem=12G -q $QUEUE -v D=$ROOT,NSPLITS=$NSPLITS_SPHERE,K=$K,c=$c,T=$T,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "SPHERE_NULL*","GROUPSC" persistEnergyGroup.sh
	done
done
