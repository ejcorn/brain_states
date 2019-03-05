#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); scan = char('$SCAN'); numClusters = [$K]; nperms = [$N]; run([basedir,'code/transprobs/nullTransitionProbs2.m']); exit"
