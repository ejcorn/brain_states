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

%% normalize states

StateMagnitude = sqrt(sum(kClusterCentroids.^2,1));
Xo = kClusterCentroids ./ StateMagnitude;		% normalize all states
Xf = Xo;	% start and end at same place, i.e. persistence energy

%% compute persistence energy for each state

DLWNullPersistenceEnergy = zeros(nperms,numClusters);

for P = 1:nperms
	disp(['Perm ',num2str(P)])
	load(fullfile(datadir,'Group_DLWNull_Gramians',['DLWNullGramianInverse',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'A_DLWnull');
	WcI_DLWnull = GRAMIAN_FAST(A_DLWnull,T);
	DLWNullPersistenceEnergy(P,:) = MIN_CONTROL_ENERGY(A_DLWnull,WcI_DLWnull,Xo,Xf,T,false); % false means don't normalize, because that is done in the previous 2 lines
end

save(fullfile(savedir,['DLWNullsPersistenceEnergy_k',num2str(numClusters),'c',num2str(c),'T',num2str(T),'.mat']),'DLWNullPersistenceEnergy');