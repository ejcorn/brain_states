#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); scan = char('$SCAN'); numClusters = [$K]; run([basedir,'code/statedynamics/nbackblockTPDur.m']); exit"
$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); scan = char('$SCAN'); numClusters = [$K]; run([basedir,'code/statedynamics/plotnbackblockTP.m']); exit"
$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); scan = char('$SCAN'); numClusters = [$K]; run([basedir,'code/statedynamics/nbackblockTPpermtest.m']); exit"