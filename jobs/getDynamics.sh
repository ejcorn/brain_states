#!/bin/bash
set -euxo pipefail

$MATPATH -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); scan = char('$SCAN'); numClusters = [$K]; run('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode/getTPDur.m'); exit"
