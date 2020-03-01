#!/bin/bash
set -euxo pipefail

cd $BD'code/control'
$RP plotEVsTP_AllTSweeps.R $D $K $BD $c
