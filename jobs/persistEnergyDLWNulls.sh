#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); nperms = [$NPERMS]; numClusters = [$K]; basedir = char('$BD'); cd(basedir); run([basedir,'code/control/persistEnergyDLWNulls.m']); exit"
