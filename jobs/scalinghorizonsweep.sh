#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); nperms = [$NPERMS]; numClusters = [$K]; nsplits = [$NSPLITS]; split = [$S]; basedir = char('$BD'); cd(basedir); run([basedir,'code/control/scalinghorizonsweep.m']); exit"