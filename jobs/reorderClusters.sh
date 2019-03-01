#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); numClusters = [$K]; run('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/kmeanscode/reorderClusters.m'); exit"
