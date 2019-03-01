#!/bin/bash
set -euxo pipefail

$MP -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); distanceMethod = char('$METHOD'); nsplits = [$N]; run([basedir,'code/assesscluster/plotzRand.m']); exit"
