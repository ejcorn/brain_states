#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); distanceMethod = char('$METHOD'); nsplits = [$N]; run('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/kmeanscode/getAssignmentspy.m'); exit"
