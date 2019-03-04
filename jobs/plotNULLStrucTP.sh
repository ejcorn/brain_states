#!/bin/bash
cd $BD'code/structp'
$RP plotSCTPNull.R $D $SCAN $THRESH $K
$RP plotNULLSubjStrucTPcorrsv2.R $D $K $THRESH
