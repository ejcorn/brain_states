#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); numClusters = [$K]; c=[$c]; basedir = char('$BD'); cd(basedir); run([basedir,'code/control/transitionEnergyDynamicsGroupv2_TSweep_CFN.m']); exit"