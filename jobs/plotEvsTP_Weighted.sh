#!/bin/bash
set -euxo pipefail

cd $BD'code/control'
$RP plotWeightedControlEnergyvsTransitionProbability.R $D $K $BD $c $NPERMS