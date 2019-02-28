#!/bin/bash
set -euxo pipefail

cd /data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode 
$RP ageDT.R $D $K
$RP ageTP.R $D $K
