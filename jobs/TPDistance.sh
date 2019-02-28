#!/bin/bash
set -euxo pipefail

cd /data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode
$RP TPdistance.R $D $K
