#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); run([basedir,'code/assesscluster/stateRepresentation.m']); exit"
$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); distanceMethod = char('$METHOD'); cd(basedir); run([basedir,'code/assesscluster/clustererror_hist.m']); exit"
