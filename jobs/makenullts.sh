#!/bin/bash
set -euxo pipefail

$MP -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); distanceMethod = char('$METHOD'); numClusters = [$K]; run([basedir,'code/assesscluster/makenullts.m']); exit"