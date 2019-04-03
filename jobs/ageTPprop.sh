#!/bin/bash
set -euxo pipefail

cd $BD'code/development'
$RP ageTPprop.R $D $K $BD
$RP ageSCTP.R $D $K $BD
