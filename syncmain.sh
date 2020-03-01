#!/bin/bash
set -euxo pipefail

rsync -avzh finalmain.sh ecornblath@chead:/data/tesla-data/ecornblath/brain_states_prepubtest/
