#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); numClusters = [$K]; run([basedir,'code/hcp/hcp_tpmat_2backnopersist.m']); exit"
