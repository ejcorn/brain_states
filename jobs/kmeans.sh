#!/bin/bash
set -euxo pipefail

$MP -nodisplay -nojvm -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); distanceMethod = char('$METHOD'); nreps = [$N]; zdim = [$Z]; split = num2str([$S]); numClusters = [$K]; run([basedir,'code/kmeans/repeatkmeans.m']); exit"
