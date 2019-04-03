#!/bin/bash
set -euxo pipefail

cd $BD'code/statedynamics/'
$RP plotRvNDur.R $D $K $BD
$RP plotRvNDwell.R $D $K $BD
$RP plotRvNRunRate.R $D $K $BD
