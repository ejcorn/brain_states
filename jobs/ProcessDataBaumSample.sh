#!/bin/bash
set -euxo pipefail
$MP -nodisplay -nojvm -r "addpath(genpath('$BD')); zdim = [$ZDIM]; scan = char('$SCAN'); lausanneScaleBOLD = [$NPARC]; extralabel = char('$LAB'); run('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/processcode/ProcessDataBaumSample.m'); exit" 
