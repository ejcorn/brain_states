#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); run([basedir,'code/control/makeGroupSC.m']); exit"
$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); run([basedir,'code/control/makeGroupSC_FA.m']); exit"