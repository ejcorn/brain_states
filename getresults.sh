#!/bin/bash
set -euxo pipefail

echo 'Enter root folder:'
#read ROOT
#ROOT='ScanCLaus250Z0reproduce_test'
ROOT='ScanCLaus250Z0prepub_test'
#ROOT='ScanCLaus250Z0final'
rsync -avzh --exclude '*SHkmeans*' ecornblath@chead:/data/tesla-data/ecornblath/brain_states_prepubtest/results/$ROOT/analyses/ results/$ROOT/
