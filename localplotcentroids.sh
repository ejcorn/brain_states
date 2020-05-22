#!/bin/bash
set -euxo pipefail
#export QT_API=pyqt
BASEDIR=~/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states
ZDIM=0
METHOD='correlation'
NPARC=250
SCAN='C'
LAB='final'

D='Scan'$SCAN'Laus'$NPARC'Z'$ZDIM$LAB


cd $BASEDIR

for K in {5..5}
do 
if [ $SCAN = 'C' ]; then
	LEN=3
else
	LEN=1
fi

for i in 3 #$(seq $LEN)  
do
	echo 'asdfsadf'
	python code/assesscluster/brainvis2.py $D $K $NPARC $SCAN $i # Fig 2a

done
python code/assesscluster/brainvis_general.py $D $K $NPARC 'centroids/motion_scrub/MotionScrubbed0.1mm' 'centroids/motion_scrub/CentroidsMotionScrubbed0.1mm_k'$K'.mat'
python code/assesscluster/brainvis_general.py $D $K $NPARC 'centroids/ipr_null/CentroidsIPRNull_k' 'centroids/ipr_null/CentroidsIPRNull_k'$K'.mat'
python code/assesscluster/brainvis_general.py $D $K $NPARC 'hcpLR/HCPCentroids_k'$K'' 'hcpLR/HCP_XHcentroids_k'$K'_R405N405'$D'.mat'
python code/assesscluster/brainvis_general.py $D $K $NPARC 'choosing_k/rvn/DownSampleTaskClusterCentroids_k' 'choosing_k/rvn/DownSampleTaskClusterCentroids_k'$K'.mat'

done

