#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); scan = char('$SCAN'); numClusters = [$K]; run('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode/nbackblockTPDur.m'); exit"

cd /data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode
$RP plotnbackblockDur.R $D $K
