# brain_states
Code to reproduce all analysis in Cornblath et al. 2018 ("Context-dependent architecture of brain state dynamics is explained by white matter connectivity and theories of network control").

Scripts are organized by their purpose in the code folder
The jobs folder contains shell scripts to allow the scripts in code folder to be submitted to a computing cluster using a Sun Grid Engine job schedule (qsub).

finalmain.sh is a bash script that will produce every figure in the paper. There are a few paths that need to be specified in this script:
