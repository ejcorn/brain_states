#!/bin/bash
set -euxo pipefail

rsync -avzh swapmain.sh ecornblath@chead:/data/tesla-data/ecornblath/brain_states/
