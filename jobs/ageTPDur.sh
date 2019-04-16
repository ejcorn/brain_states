#!/bin/bash
set -euxo pipefail

cd $BD'code/development'
$RP ageDur.R $D $K $BD
$RP ageDwell.R $D $K $BD
$RP ageBlockFO.R $D $K $BD
$RP ageBlockDwell.R $D $K $BD
$RP ageTP.R $D $K $BD
$RP ageTPintRN.R $D $K $BD
