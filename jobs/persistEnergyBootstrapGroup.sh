#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); nperms = [$NPERMS]; numClusters = [$K]; basedir = char('$BD'); cd(basedir); run([basedir,'code/control/persistEnergyBootstrapGroup.m']); exit"
cd $BD'code/control'
$RP plotEnergyDynamicsBoot.R $D $K $BD