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

Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix

Xo = kClusterCentroids(:,Xo_ind);
Xf = kClusterCentroids(:,Xf_ind);

%% compute minimum control energy for Xo-Xf

load(fullfile(datadir,['FA',num2str(lausanneScaleBOLD),'.mat']));
load(fullfile(masterdir,'analyses','control_energy',['GroupRepresentativeSC_FA_Laus',num2str(lausanneScaleBOLD),'.mat']),'G'); % load consistency threshold mask

subjectTransitionEnergy = zeros(nobs,numClusters^2);

for N = 1:nobs
	disp(['Subject ',num2str(N)])
	A = FA{N};
	A(~~eye(length(A))) = 0; % one of these matrices randomly had a NaN on diagonal
	A = A.*G; % apply consistency threshold mask
	A = NORMALIZE(A,c);
	WcI_subj = GRAMIAN_FAST(A,T); % compute gramian for whole brain control at time horizon T
	subjectTransitionEnergy(N,:) = MIN_CONTROL_ENERGY(A,WcI_subj,Xo,Xf,T,false);
end

save(fullfile(savedir,['SubjectTransitionEnergy_FA_','c',num2str(c),'T',num2str(T),'k',num2str(numClusters),'.mat']),...
	'subjectTransitionEnergy');