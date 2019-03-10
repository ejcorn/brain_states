addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy');
mkdir(savedir);

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

T = 1; %control horizon
DLWNullPersistenceEnergy = zeros(nperms,numClusters);

for P = 1:nperms
	disp(['Perm ',num2str(P)])
	load(fullfile(datadir,'Group_DLWNull_Gramians',['DLWNullGramianInverse',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'WcI_DLWnull','A_DLWnull');
	DLWNullPersistenceEnergy(P,:) = MIN_CONTROL_ENERGY(A_DLWnull,WcI_DLWnull,Xo,Xf,T,false); % false means don't normalize, see precomputeDLWNull.m for normalization
end

save(fullfile(savedir,['DLWNullsPersistenceEnergy_k',num2str(numClusters),'.mat']),'DLWNullPersistenceEnergy');