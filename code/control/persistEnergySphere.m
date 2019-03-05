addpath(genpath('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode'));
addpath(genpath('/data/tesla-data/ecornblath/matlab/brainmapping2'));
addpath(genpath('/data/tesla-data/ecornblath/matlab/BCT/'));
load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/Demographics',name_root,'.mat']);
masterdir = ['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root];
scanlab = {'RestComb','nBackComb'};

load([masterdir,'/analyses/transitionprobabilities/',scanlab{1},'TransitionProbabilities_k',num2str(numClusters),name_root,'.mat']);
load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/VolNormSC',num2str(lausanneScaleBOLD),'.mat']);
load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat']);
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;

savedir = [masterdir,'/analyses/control_energy'];
mkdir(savedir);

%% get coordinates of spatial embedding

fname = ['/data/tesla-data/ecornblath/matlab/brainmapping2/LausanneNifti/ROIv_scale',num2str(lausanneScaleBOLD),'_dilated.nii.gz'];
lausnifti = load_nii(fname);
lauscoor = zeros(nparc,3);

for i = 1:nparc
    [xind,yind,zind] = ind2sub(size(lausnifti.img),find(ismember(lausnifti.img,i)));
    lauscoor(i,:) = mean([xind,yind,zind],1);
end

D = squareform(pdist(lauscoor,'Euclidean'));

%% generate group avg SC matrix

A = zeros(nparc,nparc,nobs);
for N = 1:nobs
	A(:,:,N) = SCvolnorm{N};
end
clear SCvolnorm
         % size of connectivity matrices
A = A .* repmat(~eye(nparc),[1 1 nobs]);  % get rid of diagonal elements
frac = 1;                       % frac = 1 means your average matrix with have same # of elements as mean subject
rightnodes = 230;
hemi = [ones(rightnodes,1);2*ones(nparc - rightnodes,1)];
G = fcn_distance_dependent_threshold(A,D,hemi,frac);    % generate binary mask

% you'll need to decide how to weight G -- in the past I've done the following:
A(A == 0) = nan;    % set A = 0 to nan
W = nanmean(A,3);   % nanmean (take into account only nonzero values when computing mean edge weight)
A = G.*W;          % set G = 1 to corresponding mean edge weight
A(isnan(A)) = 0;

%% ask whether the centroid activity patterns are more stable than random activity pattern with human SC

Xo = kClusterCentroids;

%% make null centroids for each cluster that preserve spatial structure

% code from github.com/spin-test/spin-test :
% On testing for spatial correspondence between maps of human brain structure and function
% Aaron F.Alexander-Bloch et al. 2018
% https://www.sciencedirect.com/science/article/pii/S1053811918304968

% read in surface vertices of sphere for fsaverage5 and load matlab freesurfer functions
fshome = getenv('FREESURFER_HOME');
addpath(genpath(fullfile(fshome,'matlab')));
[Lvertices, Lfaces] = freesurfer_read_surf('/data/tesla-data/ecornblath/dtipreproc/freesurfer_practice/fsaverage5/surf/lh.sphere');
[Rvertices, Rfaces] = freesurfer_read_surf('/data/tesla-data/ecornblath/dtipreproc/freesurfer_practice/fsaverage5/surf/rh.sphere');
fname = ['/data/tesla-data/ecornblath/dtipreproc/freesurfer_practice/fsaverage5/label/lh.myaparc_',num2str(lausanneScaleBOLD),'.annot'];
[Lv,LL,Lct] = read_annotation(fname);
fname = ['/data/tesla-data/ecornblath/dtipreproc/freesurfer_practice/fsaverage5/label/rh.myaparc_',num2str(lausanneScaleBOLD),'.annot'];
[Rv,RL,Rct] = read_annotation(fname);    
annot.Lv = Lv; annot.LL = LL; annot.Lct = Lct; annot.Rv = Rv; annot.RL = RL; annot.Rct = Rct;
annot.Lvertices = Lvertices; annot.Lfaces = Lfaces; annot.Rvertices = Rvertices; annot.Rfaces = Rfaces;
load('/data/tesla-data/ecornblath/matlab/brainmapping2/human_regionNames.mat');
annot.roinames = roinames;

parpool(maxNumCompThreads);

nperms = 500;
Xo_Null = zeros(nparc,numClusters,nperms);
parfor K = 1:numClusters
	clusterTmp = Xo(:,K);
	[Rlabels,Llabels]= LAUS_DATA_TO_SURF_FAST(clusterTmp,annot);
	[bigrotl,bigrotr] = SpinPermuFS_EJCFAST(Llabels,Rlabels,Lvertices,Rvertices,nperms);
	%test = LAUS_SPHERE_TO_NODE(Llabels,Rlabels,nparc);		%reproduces the cluster centroid exactly
	nullTmp = LAUS_SPHERE_TO_NODE_FAST(bigrotl,bigrotr,nparc,annot);
	% leave subcortical ROIs and non-surface rendered ROIs the same as centroid
	clusterTmp = repmat(clusterTmp,[1 nperms]);
	nullTmp(nullTmp == 0) = clusterTmp(nullTmp == 0);
	nullTmp(isnan(nullTmp)) = clusterTmp(isnan(nullTmp));
	Xo_Null(:,K,:) = nullTmp;
end
cd(savedir);
save(['Xo_Null_k',num2str(numClusters),'.mat'],'Xo_Null');
disp('made null states')
% t/c replacing subcortical activity with actual subcortical activity in centroid
% maybe just find zeros and replace
StateMagnitude = sqrt(sum(kClusterCentroids.^2,1));
Xo = Xo ./ StateMagnitude;		% normalize all states
Xf = Xo;	% start and end at same place, i.e. persistence energy

NullStateMagnitude = sqrt(sum(Xo_Null.^2,1));
Xo_Null = Xo_Null ./ NullStateMagnitude;
Xf_Null = Xo_Null;


T = 1; %control horizon
Anorm = (A / max(eigs(A,1))) - eye(length(A));	% normalize A to max eig --> new max eig = 1 --> one dominant eigenmode stationary over time
% make degree preserving null model using BCT function
Arandmio = randmio_und(Anorm,10) - eye(length(A));	% for some reason randmio_und makes diagonal 0
% make degree preserving, length distribution, edge weight-length distribution preserving null model using RFB function
% Betzel & Bassett 2018 PNAS: https://doi.org/10.1073/pnas.1720186115
nbins = 11; nrewire = 1e5; [~,ArandDLW] = fcn_preserve_degseq_lengthdist(Anorm,D,nbins,nrewire); 
ArandDLW = ArandDLW + ArandDLW' - eye(length(A));


WcI = GRAMIAN(Anorm,T,false);
WcI_mio = GRAMIAN(Arandmio,T,false);
WcI_DLW = GRAMIAN(ArandDLW,T,false);
Epersist = MIN_CONTROL_ENERGY(Anorm,WcI,Xo,Xf,T,false);
Epersist_mio = MIN_CONTROL_ENERGY(Arandmio,WcI_mio,Xo,Xf,T,false);
Epersist_DLW = MIN_CONTROL_ENERGY(ArandDLW,WcI_DLW,Xo,Xf,T,false);

Epersist_Null = zeros(nperms,numClusters);
Epersist_Null_mio = zeros(nperms,numClusters);
Epersist_Null_DLW = zeros(nperms,numClusters);

for P = 1:nperms
	Epersist_Null(P,:) = MIN_CONTROL_ENERGY(Anorm,WcI,Xo_Null(:,:,P),Xf_Null(:,:,P),T,false);
	Epersist_Null_mio(P,:) = MIN_CONTROL_ENERGY(Arandmio,WcI_mio,Xo_Null(:,:,P),Xf_Null(:,:,P),T,false);
	Epersist_Null_DLW(P,:) = MIN_CONTROL_ENERGY(ArandDLW,WcI_DLW,Xo_Null(:,:,P),Xf_Null(:,:,P),T,false);
end

save(['PersistEnergySpherePerm_k_',num2str(numClusters),'.mat'],'Epersist','Epersist_Null','Epersist_mio','Epersist_Null_mio','Epersist_DLW','Epersist_Null_DLW');
