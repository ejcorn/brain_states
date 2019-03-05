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

%% calculate magnitude of states
numTrans = numClusters^2;
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numTrans); offDiag(onDiag) = [];

StateMagnitude = sqrt(sum(kClusterCentroids.^2,1));

%% calculate minimum control energy -- grp avg

tp = mean(transitionProbability,1);
Xo = kClusterCentroids ./ StateMagnitude;
Xf = kClusterCentroids ./ StateMagnitude;
%sqrt(sum(([Xo,Xf]).^2,1)) %test that all distances are same after normalization

%f = figure;		% view Xo and Xf
%subplot(1,2,1); imagesc(Xo); axis square
%subplot(1,2,2); imagesc(Xf); axis square

T = 1;		% control horizon
Anorm = (A / max(eigs(A,1))) - eye(length(A));
E = control_energy(Anorm,Xo,Xf,T,false);

[r_E_PP,p_E_PP] = corr(E',tp(onDiag)');

parpool(maxNumCompThreads - 1)
nperms = 1000;
ENull_Rewire = zeros(nperms,numClusters);
r_null_E_PP = zeros(nperms,1);
p_null_E_PP = zeros(nperms,1);
Anorm = (A / max(eigs(A,1)+1));% - eye(length(A));
parfor p = 1:nperms
	%Arand = rand(nparc); Arand = Arand + Arand';
	Arand = randmio_und(Anorm,10);% - eye(length(A));	% for some reason randmio_und makes diagonal 0
	tmp = control_energy(Arand,Xo,Xf,T,false);	% don't renormalize null matrices
	ENull_Rewire(p,:) = tmp;
	[r_null_E_PP(p),p_null_E_PP(p)] = corr(tmp',tp(onDiag)');
	disp(['Perm ',num2str(p)])
end
meanENull_Rewire = mean(ENull_Rewire,1);

rewire_pval1 = sum(E > ENull_Rewire,1) / nperms;
rewire_pval2 = sum(E < ENull_Rewire,1) / nperms;
EnergyPP_pval1 = sum(r_E_PP > r_null_E_PP) / nperms;

cd(savedir);
save(['PersistEnergyVsRewire_k_',num2str(numClusters),'.mat'], 'E','ENull_Rewire');

%% ask whether the centroid activity patterns are more stable than random activity pattern with human SC

Anorm = (A / max(eigs(A,1))) - eye(length(A));
Arand = randmio_und(Anorm,10) - eye(length(A));	% for some reason randmio_und makes diagonal 0
WcI = GRAMIAN(Anorm,T,false);
WcI_randmio = GRAMIAN(Arand,T,false);
Xo = kClusterCentroids ./ StateMagnitude;
%Xo = rand(nparc,numClusters);
Xf = Xo;
Epersist = MIN_CONTROL_ENERGY(Anorm,WcI,Xo,Xf,T,false);
Epersist_randmio = MIN_CONTROL_ENERGY(Arand,WcI_randmio,Xo,Xf,T,false);

nperms = 100;
ENull_persist = zeros(nperms,numClusters);
ENull_persist_randmio = zeros(nperms,numClusters);

tic
parfor p = 1:nperms
	Xo = kClusterCentroids(randperm(nparc),:) ./ StateMagnitude;
	Xf = Xo;
	ENull_persist(p,:) = MIN_CONTROL_ENERGY(Anorm,WcI,Xo,Xf,T,false);
	ENull_persist_randmio(p,:) = MIN_CONTROL_ENERGY(Arand,WcI_randmio,Xo,Xf,T,false);
end
toc

persist_pval1 = sum(Epersist > ENull_persist,1) / nperms;
persist_pval2 = sum(Epersist < ENull_persist,1) / nperms;
persist_randmio_pval1 = sum(Epersist_randmio > ENull_persist_randmio,1) / nperms;
persist_randmio_pval2 = sum(Epersist_randmio < ENull_persist_randmio,1) / nperms;

save(['PersistEnergyShuffleCentroid_k_',num2str(numClusters),'.mat'],'r_E_PP','r_null_E_PP','p_E_PP','p_null_E_PP','E','ENull_Rewire');

