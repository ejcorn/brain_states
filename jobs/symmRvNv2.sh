#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); scan = char('$SCAN'); numClusters = [$K]; run('$BD'code/symmRvNv2.m'); exit"

cd $BD'code/transprobs'
$RP symmRvNv2plot.R $D $K

