#!/bin/bash
set -euxo pipefail

cd $BD'code/control/'
$RP plotEnergyPersistenceProbGroup.R $D $K $BD
$RP plotEnergyDwellTimeGroup.R $D $K $BD
$RP plotPersistEnergyVsSpatialNull.R $D $K $BD