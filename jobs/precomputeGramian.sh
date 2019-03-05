#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); nsplits = [$NSPLITS]; split = [$S]; basedir = char('$BD'); cd(basedir); run([basedir,'code/control/precomputeGramian.m']); exit"
