#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); nperms = [$NPERMS]; nsplits = [$NSPLITS]; split = [$S]; basedir = char('$BD'); cd(basedir); run([basedir,'code/control/precomputeDLWNull.m']); exit"
