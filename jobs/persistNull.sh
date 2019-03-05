#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); numClusters = [$K]; run([basedir,'code/transprobs/plotPersistence.m']); exit"


