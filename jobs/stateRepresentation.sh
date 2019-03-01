#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); run([basedir,'code/assesscluster/stateRepresentation.m']); exit"
