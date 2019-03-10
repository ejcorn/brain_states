addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy');
mkdir(savedir);

%% load states

% empirical states
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;

%% normalize states

StateMagnitude = sqrt(sum(kClusterCentroids.^2,1));
Xo = kClusterCentroids ./ StateMagnitude;		% normalize all states
Xf = Xo;	% start and end at same place, i.e. persistence energy

%% compute minimum control energy for Xo-Xf

load(fullfile(datadir,['VolNormSC',num2str(lausanneScaleBOLD),'.mat']));
T = 1; %control horizon

subjectPersistenceEnergy = zeros(nobs,numClusters);

for N = 1:nobs
	disp(['Subject ',num2str(N)])
	load(fullfile(datadir,'Gramians',['GramianInverse',num2str(lausanneScaleBOLD),'Subject',num2str(N),'.mat']),'WcI_subj');
	A = SCvolnorm{N};
	A(~~eye(length(A))) = 0; % one of these matrices randomly had a NaN on diagonal
	subjectPersistenceEnergy(N,:) = MIN_CONTROL_ENERGY(A,WcI_subj,Xo,Xf,T,true);
end

save(fullfile(savedir,['SubjectPersistenceEnergy_k',num2str(numClusters),'.mat']),'subjectPersistenceEnergy');