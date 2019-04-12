#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); nsplits = [$NSPLITS]; numClusters = [$K]; c=[$c]; T=[$T]; run([basedir,'code/control/persistEnergyGroup.m']); exit"
cd $BD'code/control/'
$RP plotEnergySphereGroup.R $D $K $BD $c $T