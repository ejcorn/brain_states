#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); numClusters = [$K]; c=[$c]; InputSystems=[$SYS]+1; SystemLabel = char('$SYSNAME'); basedir = char('$BD'); cd(basedir); run([basedir,'code/control/transitionEnergyPredictedTaskInputs_CFN.m']); exit"