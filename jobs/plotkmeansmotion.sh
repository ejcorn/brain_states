#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); distanceMethod = char('$METHOD'); numClusters = [$K]; run([basedir,'code/assesscluster/motioncluster_plot.m']); exit"