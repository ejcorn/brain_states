# brain_states
Code to reproduce all analysis in Cornblath et al. 2020 ("Temporal sequences of brain activity at rest are constrained by white matter structure and modulated by cognitive demands").

The most simple way to re-implement this work is to start by looking at the script `example.m` in the master folder before digging into the full set of scripts called by `finalmain.sh`.

## Requirements:
  - MATLAB R2017a or later
  - R 3.2.5 or later, with packages:
    - ggplot2
    - R.matlab
    - RColorBrewer
    - lm.beta
    - reshape2
    - viridis
    - plotrix
  - Hardware: a computing cluster using Sun Grid Engine job scheduler, ability to request cores with at least 16G of RAM
  
This software was tested on GNU Linux using the Center for Functional Neuroimaging computing cluster (https://cfn.upenn.edu/).

## Directory structure and path specification

Scripts are organized by their purpose in the code folder. The jobs folder contains shell scripts to allow the scripts in code folder to be submitted to a computing cluster using a Sun Grid Engine job scheduler (qsub). 

finalmain.sh is a bash script that will produce every figure in the paper. There are a few paths that need to be specified in this script:
  - BASEDIR: the path to the master branch of this repository, i.e. /Users/Eli/Dropbox/brain_states
  - MATPATH: the path to the user's MATLAB binary, i.e. /Applications/MATLAB_R2017a.app/bin/matlab
  - RPATH: the path to the user's Rscript function i.e. Rscript

Additionally, data dependencies are as follows (where `$NPARC` is the parcellation scale and `$ID` is a subject ID):
  - `datadir_main`: at the top of processDataBaumSample.m, specifies the path to processed data. Within this directory there should be:
      - `subject_demographics/*.csv` with csv file of demographics, containing scanids and bblids
      - `diffusion_data/volNormSC/Lausannne${NPARC}/${ID}_volNormStreamline_LausanneScale${NPARC}.mat`
      - `diffusion_data/FA/Lausanne${NPARC}/${ID}_FA_LausanneScale${NPARC}.mat`
      - `rest_data/Lausanne${NPARC}/${ID}_Lausanne${NPARC}_ts.1D`
      - `nback_data/Lausanne${NPARC}/${ID}_Lausanne${NPARC}_ts.1D`
  - `data/*`: the main data folder of the pipeline
      - .csv files with n-back performance data
      - .csv with n-back framewise motion time series
      - `motion_files/*`: contains `\n` delimited text files named `${BBLID}_${SCANID}_relRMS.1D` with time series of framewise motion for each subject's resting state scan
  
ProcessDataBaumSample.m requires file paths to the resting state and n-back BOLD time series, as well as structural adjacency matrices from diffusion tractography. To demo this code without obtaining the necessary BOLD data, one could replace the variables "concTS" and "SCVolnorm" in this script with random numbers or your own BOLD data.

code/miscfxns/addpaths.m also requires specification of path to the BCT:
  - Brain Connectivity Toolbox, MATLAB version: https://sites.google.com/site/bctnet/

## Input specification

The user can also specify certain parameters in finalmain.sh:

  - ZDIM: should time series be z-scored? We did not z-score but one can if desired.
  - METHOD: specify a valid distance function for MATLAB kmeans. We used 'correlation'.
  - NPARC: specify Lausanne parcellation scale. Integer, either 60, 125, or 250.
  - SCAN: 'C' uses both rest and n-back. 'R' and 'N' use rest or n-back only, but this functionality is deprecated
  - LAB: string, any extra label you'd like to place on the output
  - NSPLITS: dictates how many repetitions of k-means clustering to perform for each value of K. For the submitted manuscript we used NSPLITS=100. Empirically, you will safely arrive at the solution presented in the main text for K=5 using NSPLITS=20.
  - K: number of clusters for which all analyses will be performed.
 
finalmain.sh will submit NSPLITS-times-10 + ~600 jobs. In our experience, this will take about a day to run depending on the availability of computing cores. We have not attempted to run this on a normal desktop computer.

## Understanding individual scripts in context

Most of the main analysis scripts are meant to load output from upstream scripts. The file finalmain.sh calls all individual scripts to a job scheduler that will run each script in order of dependency. The initial part of the pipeline collects processed BOLD data and subject information, performs clustering, and saves the cluster assignments. Three variables are passed into almost every MATLAB script:

  -`basedir`: home/working directory for the project see (**Directory structure and path specification**)
  -`name_root`: name of specific output folder for given set **Input specification**
  -`numClusters`: number of clusters for clustering solution you want to analyze

Almost every MATLAB script begins by loading general data (parcellation, number of subjects) and clustering related data (cluster assignments, cluster centroids, cluster names, assignment of TRs to scan or subject). Some scripts will also load the .csv file with BOLD data. So when reading the scripts, expect that some loaded variables are in the workspace but are only defined in other scripts.

In order to run just an individual script without rerunning everything, you can paste into your console predefined values for the 3 variables above and it should be sufficient to run most scripts provided the dependent data has already been generated. An example of these "start up" commands is in the file startup_params.txt

Please contact Eli Cornblath (Eli ~`DOT`~ Cornblath ~`AT`~ pennmedicine.upenn.edu) with any questions regarding this code.
