#!/bin/bash
set -euxo pipefail

cd $BD'code/control/'
$RP plotEnergySphere.R $D $K $BD
#$RP plotEnergyDynamics.R $D $K $BD
