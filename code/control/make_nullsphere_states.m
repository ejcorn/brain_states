addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy',['nullstates_k',num2str(numClusters)]);
mkdir(savedir);

% make sure you get different states each iteration
rng(str2num([num2str(numClusters),num2str(split_i)]));

%% make null centroids for each cluster that preserve spatial structure

% code from github.com/spin-test/spin-test :
% On testing for spatial correspondence between maps of human brain structure and function
% Aaron F.Alexander-Bloch et al. 2018
% https://www.sciencedirect.com/science/article/pii/S1053811918304968

% read in surface vertices of sphere for fsaverage5 and load matlab freesurfer functions

load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;

% get all information from surface files into struct, which functions below will use to carry out sphere-based permutation

fshome = getenv('FREESURFER_HOME');
addpath(genpath(fullfile(fshome,'matlab')));
[Lvertices, Lfaces] = freesurfer_read_surf(fullfile(basedir,'data/annot/lh.sphere'));
[Rvertices, Rfaces] = freesurfer_read_surf(fullfile(basedir,'data/annot/lh.sphere'));
fname = fullfile(basedir,['data/annot/lh.myaparc_',num2str(lausanneScaleBOLD),'.annot']);
[Lv,LL,Lct] = read_annotation(fname);
fname = fullfile(basedir,['data/annot/rh.myaparc_',num2str(lausanneScaleBOLD),'.annot']);
[Rv,RL,Rct] = read_annotation(fname);    
annot.Lv = Lv; annot.LL = LL; annot.Lct = Lct; annot.Rv = Rv; annot.RL = RL; annot.Rct = Rct;
annot.Lvertices = Lvertices; annot.Lfaces = Lfaces; annot.Rvertices = Rvertices; annot.Rfaces = Rfaces;
load('data/human_regionNames.mat');
annot.roinames = roinames;


%cluster_i is a variable passed in through bash script to parallelize this computation by computing one cluster at a time

nperms = 20;

Xo_Null = zeros(nparc,nperms);		% make matrix to hold permuted centroid

clusterTmp = kClusterCentroids(:,cluster_i);
disp('converting node data to surface')
[Rlabels,Llabels]= LAUS_DATA_TO_SURF_FAST(clusterTmp,annot);
disp(['rotating sphere ',num2str(nperms),' times'])
[bigrotl,bigrotr] = SpinPermuFS_EJCFAST(Llabels,Rlabels,Lvertices,Rvertices,nperms);
%test = LAUS_SPHERE_TO_NODE(Llabels,Rlabels,nparc);		%reproduces the cluster centroid exactly
disp(['applying rotation'])
nullTmp = LAUS_SPHERE_TO_NODE_FAST(bigrotl,bigrotr,nparc,annot);
% leave subcortical ROIs and non-surface rendered ROIs the same as centroid
clusterTmp = repmat(clusterTmp,[1 nperms]);
nullTmp(nullTmp == 0) = clusterTmp(nullTmp == 0);	% replace subcortical ROIs (0's as output of LAUS_SPHERE_TO_NODE_FAST) with existing activity levels
nullTmp(isnan(nullTmp)) = clusterTmp(isnan(nullTmp));
Xo_Null = nullTmp;
save(fullfile(savedir,['Xo_Null_Cluster',num2str(cluster_i),'_k',num2str(numClusters),'Split',num2str(split_i),'.mat']),'Xo_Null');
disp('saved null states')

% t/c replacing subcortical activity with actual subcortical activity in centroid
% maybe just find zeros and replace

