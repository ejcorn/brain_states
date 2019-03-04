#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); numClusters = [$K]; run([basedir,'code/hcp/hcp_cluster.m']); exit"