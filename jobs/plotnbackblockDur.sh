#!/bin/bash
set -euxo pipefail

cd $BD'code/statedynamics'
$RP plotnbackblockDur.R $D $K $BD
$RP plotnbackblockDwell.R $D $K $BD
