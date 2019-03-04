#!/bin/bash
set -euxo pipefail

cd $BD'code/transprobs'
$RP TPdistance.R $D $K $BD
