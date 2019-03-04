#!/bin/bash
set -euxo pipefail

cd $BD'code/development'
$RP ageTPprop.R $D $K $BD
