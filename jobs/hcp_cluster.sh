#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); distanceMethod = char('$METHOD'); basedir = char('$BD'); cd(basedir); numClusters = [$K]; run([basedir,'code/hcp/hcp_cluster.m']); exit"
