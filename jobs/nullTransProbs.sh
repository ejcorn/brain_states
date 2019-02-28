#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); scan = char('$SCAN'); numClusters = [$K]; nperms = [$N]; run('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode/nullTransitionProbs2.m'); exit"
