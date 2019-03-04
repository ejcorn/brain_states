#!/bin/bash
set -euxo pipefail

cd $BD'code/development'
$RP ageDur.R $D $K $BD
$RP ageTP.R $D $K $BD
