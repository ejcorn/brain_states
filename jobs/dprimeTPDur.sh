#!/bin/bash
set -euxo pipefail

cd $BD'code/cognition'
$RP dprimeTPDur.R $D $K $BD
$RP dprimeTPscatter.R $D $K $BD
$RP dprimeDwell.R $D $K $BD
