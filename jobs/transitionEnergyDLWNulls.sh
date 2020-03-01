#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); nperms = [$NPERMS]; numClusters = [$K]; c=[$c]; T=[$T]; basedir = char('$BD'); cd(basedir); run([basedir,'code/control/transitionEnergyDLWNulls.m']); exit"
$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); nperms = [$NPERMS]; numClusters = [$K]; c=[$c]; T=[$T]; basedir = char('$BD'); cd(basedir); run([basedir,'code/control/plot_transitionEnergyMagnitudeVsDLWNull.m']); exit"