#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); numClusters = [$K]; thrsh = [$THRESH]; run([basedir,'code/structp/SCSTPtrans.m']); exit"
