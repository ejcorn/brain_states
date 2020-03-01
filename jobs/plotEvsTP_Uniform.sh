#!/bin/bash
set -euxo pipefail

cd $BD'code/control'
$RP plotEnergyVsTransitionProbability.R $D $K $BD $c