#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); lausanneScaleBOLD = [$NPARC]; scan = char('$SCAN'); run([basedir,'code/process/rickNullSC.m']); exit"

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); lausanneScaleBOLD = [$NPARC]; scan = char('$SCAN'); run([basedir,'code/process/BCTSCnull.m']); exit"
