# brain_states
Code to reproduce all analysis in Cornblath et al. 2018 ("Context-dependent architecture of brain state dynamics is explained by white matter connectivity and theories of network control").

Requirements:
  - MATLAB R2017a or later
  - R 3.2.5 or later, with packages:
    - ggplot2, R.matlab, RColorBrewer

Scripts are organized by their purpose in the code folder. The jobs folder contains shell scripts to allow the scripts in code folder to be submitted to a computing cluster using a Sun Grid Engine job schedule (qsub).

finalmain.sh is a bash script that will produce every figure in the paper. There are a few paths that need to be specified in this script:

  #BASEDIR: the path to the master branch of this repository, i.e. /Users/Eli/Dropbox/brain_states
  #MATPATH: the path to the user's MATLAB binary, i.e. /Applications/MATLAB_R2017a.app/bin/matlab
  #RPATH: the path to the user's Rscript function i.e. Rscript

The user can also specify certain parameters in finalmain.sh:

  -ZDIM: should time series be z-scored? We did not z-score but one can if desired.
  -METHOD: specify a valid distance function for MATLAB kmeans. We used 'correlation'.
  -NPARC: specify Lausanne parcellation scale. Integer, either 60, 125, or 250.
  -SCAN: 'C' uses both rest and n-back. 'R' and 'N' use rest or n-back only, but this functionality is deprecated
  -LAB: string, any extra label you'd like to place on the output
  
