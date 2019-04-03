#!/bin/bash
set -euxo pipefail

cd $BD'code/structp'
$RP plotSubjStrucTPcorrsv2.R $D $K $THRESH $BD
