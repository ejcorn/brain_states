#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); nperms = [$NPERMS]; nsplits = [$NSPLITS]; split = [$S]; numClusters = [$K]; c=[$c]; T=[$T]; SystemLabel = char('$SYSTEM'); ControlLabel = char('$CONTROL'); basedir = char('$BD'); cd(basedir); run([basedir,'code/control/transitionEnergySystemWeightedControlBCTNull.m']); exit"