#!/bin/bash
set -euxo pipefail

echo 'Enter root folder:'
read ROOT
rsync -avzh ecornblath@chead:/data/tesla-data/ecornblath/brain_states/results/$ROOT/analyses/ results/$ROOT/