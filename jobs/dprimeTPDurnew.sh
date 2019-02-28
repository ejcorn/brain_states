#!/bin/bash
set -euxo pipefail

cd /data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode
$RP dprimeTPDur.R $D $K
$RP dprimeTPscatter.R $D $K

