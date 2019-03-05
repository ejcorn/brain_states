#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); cluster_i = [$Ki]; numClusters = [$K]; run([basedir,'code/control/make_nullsphere_states.m']); exit"