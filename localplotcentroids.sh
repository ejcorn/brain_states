#!/bin/bash
set -euxo pipefail
#export QT_API=pyqt
BASEDIR=~/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states
ZDIM=0
METHOD='correlation'
NPARC=250
SCAN='C'
LAB='corrfinal'

D='Scan'$SCAN'Laus'$NPARC'Z'$ZDIM$LAB


cd $BASEDIR

for K in {5..5}
do 
if [ $SCAN = 'C' ]; then
	LEN=3
else
	LEN=1
fi

for i in 3 
#$(seq $LEN)  
do

	python code/assesscluster/brainvis2.py $D $K $NPARC $SCAN $i $BASEDIR

done
done
