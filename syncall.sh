#!/bin/bash
set -euxo pipefail

rsync -avzh . --exclude 'data' --exclude 'all_*' --exclude 'results' ecornblath@chead:/data/tesla-data/ecornblath/brain_states_prepubtest/