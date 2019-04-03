#!/bin/bash
set -euxo pipefail

cd $BD'code/structp'
$RP plotSCTP.R $D $SCAN $THRESH $K $BD
