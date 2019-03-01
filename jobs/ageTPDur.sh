#!/bin/bash
set -euxo pipefail

cd $BD'code/development'
$RP ageDur.R $D $K
$RP ageTP.R $D $K
