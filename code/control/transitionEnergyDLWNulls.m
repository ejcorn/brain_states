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

% empirical states
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

%% define state transitions

Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix

Xo = kClusterCentroids(:,Xo_ind);
Xf = kClusterCentroids(:,Xf_ind);

%% compute transition energy for each state

DLWNullTransitionEnergy = zeros(nperms,numClusters^2);

for P = 1:nperms
	disp(['Perm ',num2str(P)])
	load(fullfile(datadir,'Group_DLWNull_Gramians',['DLWNull_FA_Laus',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'A_DLWnull');
	A_DLWnullnorm = NORMALIZE(A_DLWnull,c); % stabilize null matrix
	WcI_DLWnull = GRAMIAN_FAST(A_DLWnullnorm,T);
	DLWNullTransitionEnergy(P,:) = MIN_CONTROL_ENERGY(A_DLWnullnorm,WcI_DLWnull,Xo,Xf,T,false); % false means don't normalize
end

save(fullfile(savedir,['DLWNullsTransitionEnergy_FA_k',num2str(numClusters),'c',num2str(c),'T',num2str(T),'.mat']),'DLWNullTransitionEnergy');