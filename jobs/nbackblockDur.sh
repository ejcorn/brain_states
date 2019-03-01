#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); scan = char('$SCAN'); numClusters = [$K]; run([basedir,'code/statedynamics/nbackblockTPDur.m']); exit"

cd /data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode
$RP plotnbackblockDur.R $D $K
