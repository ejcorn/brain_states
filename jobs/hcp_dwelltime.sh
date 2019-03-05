#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); numClusters = [$K]; run([basedir,'code/hcp/hcp_dwelltime.m']); exit"

cd $BD'code/hcp'
$RP HCPvsPNCDwell.R $D $K $BD