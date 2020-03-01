addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy');
mkdir(savedir);

% variables that get passed in:
% c: normalization factor
% normalize each matrix by 1/c+max(eig(A)), then subtract identity. ensures max eigenvalue is [0,-1), such that matrices are stable and 
% at long time horizons all nodes either go to 0 or to a constant non-zero eigenvector associated with max eig
% T: control horizon, time over which to exert control

%% load states

% null states
load(fullfile(savedir,['nullstates_k',num2str(numClusters)],['Xo_Null_Cluster1_k',num2str(numClusters),'Split1.mat']),'Xo_Null');
nreps = size(Xo_Null,2);
nperms = nreps*nsplits;

Xo_Null_all = zeros(nparc,numClusters,nperms);	% preallocate 3D array to store every null state for EVERY cluster

for K = 1:numClusters
	Xo_Null_k = zeros(nparc,nperms);	% preallocate 2D array to store every null state for EACH cluster
	for S = 1:nsplits
		load(fullfile(savedir,['nullstates_k',num2str(numClusters)],['Xo_Null_Cluster',num2str(K),'_k',num2str(numClusters),'Split',num2str(S),'.mat']),'Xo_Null');
		Xo_Null_k(:,(1+nreps*(S-1)):(nreps*S)) = Xo_Null;	% concatenate splits into one matrix
	end
	Xo_Null_all(:,K,:) = Xo_Null_k;	% concatenate null states from all clusters into one matrix
end

% empirical states
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
%% normalize states

StateMagnitude = sqrt(sum(kClusterCentroids.^2,1));
Xo = kClusterCentroids ./ StateMagnitude;		% normalize all states
Xf = Xo;	% start and end at same place, i.e. persistence energy

NullStateMagnitude = sqrt(sum(Xo_Null_all.^2,1));
Xo_Null_all = Xo_Null_all ./ NullStateMagnitude;	% normalize null states
Xf_Null = Xo_Null_all;

%% load structure, normalize, generate nulls
load(fullfile(savedir,['GroupRepresentativeSC_FA_Laus',num2str(lausanneScaleBOLD),'.mat']),'A','D');
%load(fullfile(datadir,'BootstrapGroupSC',['BootstrapGroupRepresentativeSC',num2str(lausanneScaleBOLD),'Perm1.mat']),'WcI_boot','G_boot','A_boot','bootsamp');

A = NORMALIZE(A,c);
% make degree preserving null model using BCT function
Arandmio = randmio_und(A+eye(length(A)),10);	% need to make diagonal 0 for randmio_und -- doesn't deal w negatives
Arandmio = Arandmio - eye(length(A)); % subtract identity back out to ensure lambda max (-1,0]

% make degree preserving, length distribution, edge weight-length distribution preserving null model using RFB function
% Betzel & Bassett 2018 PNAS: https://doi.org/10.1073/pnas.1720186115
nbins = 11; nrewire = 1e5; [~,ArandDLW] = fcn_preserve_degseq_lengthdist(A,D,nbins,nrewire); 
% this function outputs an asymmetric matrix with half of the global network strength as input matrix
% by global network strength I mean sum of all connections
% it also will set diagonal to 0 by default
% need to symmetrize and double global strength, then subtract identity in order to compare
ArandDLW = ArandDLW + ArandDLW' - eye(length(A));	

%% compute minimum control energy required to maintain each state

% compute inverse of controllability gramian for each network

WcI = GRAMIAN_FAST(A,T);
WcI_mio = GRAMIAN_FAST(Arandmio,T);
WcI_DLW = GRAMIAN_FAST(ArandDLW,T);

% compute minimum control energy required to maintain cluster centers
Epersist = MIN_CONTROL_ENERGY(A,WcI,Xo,Xf,T,false);
Epersist_mio = MIN_CONTROL_ENERGY(Arandmio,WcI_mio,Xo,Xf,T,false);
Epersist_DLW = MIN_CONTROL_ENERGY(ArandDLW,WcI_DLW,Xo,Xf,T,false);

% compute minimum control energy required to maintain sphere-permuted null cluster centers
Epersist_Null = zeros(nperms,numClusters);
Epersist_Null_mio = zeros(nperms,numClusters);
Epersist_Null_DLW = zeros(nperms,numClusters);

for P = 1:nperms
	disp(['Null state ',num2str(P)])
	Epersist_Null(P,:) = MIN_CONTROL_ENERGY(A,WcI,Xo_Null_all(:,:,P),Xf_Null(:,:,P),T,false);
	Epersist_Null_mio(P,:) = MIN_CONTROL_ENERGY(Arandmio,WcI_mio,Xo_Null_all(:,:,P),Xf_Null(:,:,P),T,false);
	Epersist_Null_DLW(P,:) = MIN_CONTROL_ENERGY(ArandDLW,WcI_DLW,Xo_Null_all(:,:,P),Xf_Null(:,:,P),T,false);
end

save(fullfile(savedir,['PersistEnergySpherePerm_FA_k',num2str(numClusters),'c',num2str(c),'T',num2str(T),'.mat']),'Epersist','Epersist_Null','Epersist_mio','Epersist_Null_mio','Epersist_DLW','Epersist_Null_DLW');