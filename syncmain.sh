#!/bin/bash
set -euxo pipefail

rsync -avzh testmain2.sh ecornblath@chead:/data/tesla-data/ecornblath/brain_states/
