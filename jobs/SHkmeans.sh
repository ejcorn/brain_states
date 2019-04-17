#!/bin/bash
set -euxo pipefail

$MP -nodisplay -nojvm -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); split = [$S]; cd(basedir); distanceMethod = char('$METHOD'); zdim=[$Z]; numClusters = [$K]; run([basedir,'code/kmeans/SHrepeatkmeans.m']); exit"