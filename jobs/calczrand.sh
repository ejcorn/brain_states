#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); numClusters = [$K]; distanceMethod = char('$METHOD'); nsplits = [$N]; run([basedir,'code/assesscluster/calczRand.m']); exit"
