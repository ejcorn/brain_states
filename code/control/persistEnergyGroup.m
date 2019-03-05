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
		load(fullfile(savedir,['Xo_Null_Cluster',num2str(K),'_k',num2str(numClusters),'Split',num2str(S),'.mat']),'Xo_Null');
		Xo_Null_k(:,(1+nreps*(S-1)):(nreps*S)) = Xo_Null;	% concatenate splits into one matrix
	end
	Xo_Null_all(:,K,:) = Xo_Null_k;	% concatenate null states from all clusters into one matrix
end

% empirical states
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;

%% normalize states

StateMagnitude = sqrt(sum(kClusterCentroids.^2,1));
Xo = Xo ./ StateMagnitude;		% normalize all states
Xf = Xo;	% start and end at same place, i.e. persistence energy

NullStateMagnitude = sqrt(sum(Xo_Null_all.^2,1));
Xo_Null_all = Xo_Null_all ./ NullStateMagnitude;	% normalize null states
Xf_Null = Xo_Null_all;

%% load structure, normalize, generate nulls
load(fullfile(savedir,['GroupRepresentativeSC_Laus',num2str(lausanneScaleBOLD),'.mat']),'A');

T = 1; %control horizon
Anorm = (A / max(eigs(A,1))) - eye(length(A));	% normalize A to max eig --> new max eig = 1 --> one dominant eigenmode stationary over time
% make degree preserving null model using BCT function
Arandmio = randmio_und(Anorm,10) - eye(length(A));	% for some reason randmio_und makes diagonal 0
% make degree preserving, length distribution, edge weight-length distribution preserving null model using RFB function
% Betzel & Bassett 2018 PNAS: https://doi.org/10.1073/pnas.1720186115
nbins = 11; nrewire = 1e5; [~,ArandDLW] = fcn_preserve_degseq_lengthdist(Anorm,D,nbins,nrewire); 
ArandDLW = ArandDLW + ArandDLW' - eye(length(A));

%% compute minimum control energy required to maintain each state

% compute inverse of controllability gramian for each network
WcI = GRAMIAN(Anorm,T,false);
WcI_mio = GRAMIAN(Arandmio,T,false);
WcI_DLW = GRAMIAN(ArandDLW,T,false);

% compute minimum control energy required to maintain cluster centers
Epersist = MIN_CONTROL_ENERGY(Anorm,WcI,Xo,Xf,T,false);
Epersist_mio = MIN_CONTROL_ENERGY(Arandmio,WcI_mio,Xo,Xf,T,false);
Epersist_DLW = MIN_CONTROL_ENERGY(ArandDLW,WcI_DLW,Xo,Xf,T,false);

% compute minimum control energy required to maintain sphere-permuted null cluster centers
Epersist_Null = zeros(nperms,numClusters);
Epersist_Null_mio = zeros(nperms,numClusters);
Epersist_Null_DLW = zeros(nperms,numClusters);

for P = 1:nperms
	Epersist_Null(P,:) = MIN_CONTROL_ENERGY(Anorm,WcI,Xo_Null_all(:,:,P),Xf_Null(:,:,P),T,false);
	Epersist_Null_mio(P,:) = MIN_CONTROL_ENERGY(Arandmio,WcI_mio,Xo_Null_all(:,:,P),Xf_Null(:,:,P),T,false);
	Epersist_Null_DLW(P,:) = MIN_CONTROL_ENERGY(ArandDLW,WcI_DLW,Xo_Null_all(:,:,P),Xf_Null(:,:,P),T,false);
end

save(fullfile(savedir,['PersistEnergySpherePerm_k_',num2str(numClusters),'.mat']),'Epersist','Epersist_Null','Epersist_mio','Epersist_Null_mio','Epersist_DLW','Epersist_Null_DLW');
