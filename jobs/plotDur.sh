#!/bin/bash
set -euxo pipefail

cd $BD'code/statedynamics/'
$RP plotRvNDur.R $D $K $
