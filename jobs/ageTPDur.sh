#!/bin/bash
set -euxo pipefail

cd $BD'code/development'
$RP ageDwellRestBlocks.R $D $K $BD
$RP ageFORestBlocks.R $D $K $BD
$RP ageBlockTP.R $D $K $BD
