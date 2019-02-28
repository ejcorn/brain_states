#!/bin/bash
cd /data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode
$RP plotSCTPNull.R $D $SCAN $THRESH $K
$RP plotNULLSubjStrucTPcorrsv2.R $D $K $THRESH
