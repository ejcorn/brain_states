#!/bin/bash
set -euxo pipefail
$MP -nodisplay -nojvm -r "addpath(genpath('$BD')); zdim = [$ZDIM]; scan = char('$SCAN'); lausanneScaleBOLD = [$NPARC]; extralabel = char('$LAB'); run('$BD/code/process/ProcessDataBaumSample.m'); exit" 
