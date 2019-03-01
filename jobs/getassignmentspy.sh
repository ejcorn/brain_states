#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); distanceMethod = char('$METHOD'); nsplits = [$N]; run([basedir,'code/kmeans/getAssignmentspy.m']); exit"
