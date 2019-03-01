#!/bin/bash
set -euxo pipefail

$MP -nodisplay -nojvm -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); distanceMethod = char('$METHOD'); nreps = [$N]; zdim = [$Z]; split = num2str([$S]); numClusters = [$K]; run('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/kmeanscode/repeatkmeans.m'); exit"
