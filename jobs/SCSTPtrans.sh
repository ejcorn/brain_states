#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); lausanneScaleBOLD = [$NPARC]; scan = char('$SCAN'); numClusters = [$K]; thrsh = [$THRESH]; run('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode/SCSTPtrans.m'); exit"
