#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); numClusters = [$K]; c=[$c]; T=[$T]; run([basedir,'code/control/transitionEnergySubjectsv2_FA.m']); exit"

cd $BD'code/development'
$RP ageTransitionEnergy.R $D $K $BD $c $T