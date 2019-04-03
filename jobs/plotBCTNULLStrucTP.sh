#!/bin/bash
set -euxo pipefail

cd $BD'code/structp'
$RP plotBCTNullSubjStrucTPcorrsv2.R $D $K $THRESH $BD
