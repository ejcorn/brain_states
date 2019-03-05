#!/bin/bash
set -euxo pipefail
$MP -nodisplay -nojvm -r "addpath(genpath('$BD')); basedir = char('$BD'); cd(basedir); zdim = [$ZDIM]; scan = char('$SCAN'); lausanneScaleBOLD = [$NPARC]; extralabel = char('$LAB'); run([basedir,'code/process/ProcessDataBaumSample.m']); exit" 
