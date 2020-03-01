#!/bin/bash
set -euxo pipefail

cd $BD'code/cognition'
$RP dprimeTPbyBlock.R $D $K $BD
$RP dprimeBlockFO.R $D $K $BD
$RP dprimeDwell.R $D $K $BD

$MP -nodisplay -r "addpath(genpath('$BD')); name_root = char('$D'); basedir = char('$BD'); cd(basedir); numClusters = [$K]; run([basedir,'code/cognition/dprimeTP2Back.m']); exit"