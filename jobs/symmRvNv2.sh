#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); scan = char('$SCAN'); numClusters = [$K]; run([basedir,'code/transprobs/symmRvNv2.m']); exit"

cd $BD'code/transprobs'
$RP symmRvNv2plot.R $D $K $BD

