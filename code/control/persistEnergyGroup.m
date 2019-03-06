addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy');
mkdir(savedir);
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
load(fullfile(savedir,['GroupRepresentativeSC_Laus',num2str(lausanneScaleBOLD),'.mat']),'A','D');
c = 1;
Anorm = (A / max(eig(A))) - c*eye(length(A));   % normalize by max eigenvalue
max(eig(Anorm)) % max eigenvalue should be v. close to 1 - c;

T = 1; %control horizon
%Anorm = (A / max(eigs(A,1))) - eye(length(A));	% normalize A to max eig --> new max eig = 0 --> one dominant eigenmode stationary over time
% make degree preserving null model using BCT function
Arandmio = randmio_und(A,10);	% for some reason randmio_und makes diagonal 0
% make degree preserving, length distribution, edge weight-length distribution preserving null model using RFB function
% Betzel & Bassett 2018 PNAS: https://doi.org/10.1073/pnas.1720186115
nbins = 11; nrewire = 1e5; [~,ArandDLW] = fcn_preserve_degseq_lengthdist(A,D,nbins,nrewire); 
ArandDLW = ArandDLW + ArandDLW';

%% compute minimum control energy required to maintain each state

% compute inverse of controllability gramian for each network
% normalize all networks to max eig then subtract identity to ensure max eigval = 0 so system converges on single eigenmode over time
WcI = GRAMIAN(A,T,true);
WcI_mio = GRAMIAN(Arandmio,T,true);
WcI_DLW = GRAMIAN(ArandDLW,T,true);

% compute minimum control energy required to maintain cluster centers
Epersist = MIN_CONTROL_ENERGY(A,WcI,Xo,Xf,T,true);
Epersist_mio = MIN_CONTROL_ENERGY(Arandmio,WcI_mio,Xo,Xf,T,true);
Epersist_DLW = MIN_CONTROL_ENERGY(ArandDLW,WcI_DLW,Xo,Xf,T,true);

% compute minimum control energy required to maintain sphere-permuted null cluster centers
Epersist_Null = zeros(nperms,numClusters);
Epersist_Null_mio = zeros(nperms,numClusters);
Epersist_Null_DLW = zeros(nperms,numClusters);

for P = 1:nperms
	disp(['Null state ',num2str(P)])
	Epersist_Null(P,:) = MIN_CONTROL_ENERGY(A,WcI,Xo_Null_all(:,:,P),Xf_Null(:,:,P),T,true);
	Epersist_Null_mio(P,:) = MIN_CONTROL_ENERGY(Arandmio,WcI_mio,Xo_Null_all(:,:,P),Xf_Null(:,:,P),T,true);
	Epersist_Null_DLW(P,:) = MIN_CONTROL_ENERGY(ArandDLW,WcI_DLW,Xo_Null_all(:,:,P),Xf_Null(:,:,P),T,true);
end

save(fullfile(savedir,['PersistEnergySpherePerm_k_',num2str(numClusters),'.mat']),'Epersist','Epersist_Null','Epersist_mio','Epersist_Null_mio','Epersist_DLW','Epersist_Null_DLW');
