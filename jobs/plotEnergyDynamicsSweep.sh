#!/bin/bash
set -euxo pipefail

cd $BD'code/control'
$RP plotEnergyDynamicsSweep.R $D $K $BD