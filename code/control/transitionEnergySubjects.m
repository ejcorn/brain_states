addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy');
mkdir(savedir);

%% load states

% empirical states
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;

%% define sets of initial and final conditions
% for now, not normalizing since I'm using it for comparisons across people not between states

i = repmat(1:numClusters,[1 numClusters]);	% define initial states
Xo = kClusterCentroids(:,i);
j = repelem(1:numClusters,numClusters);		% define final states
Xf = kClusterCentroids(:,j);
% order of Xo and Xf is equiv. to flattening by row, like in GET_TRANS_PROBS.m

sqrt(sum((Xo-Xf).^2,1)) %test all distances between states
%% compute minimum control energy for Xo-Xf

load(fullfile(datadir,['VolNormSC',num2str(lausanneScaleBOLD),'.mat']));
T = 1; %control horizon

subjectTransitionEnergy = zeros(nobs,numClusters^2);

for N = 1:nobs
	disp(['Subject ',num2str(N)])
	load(fullfile(datadir,'Gramians',['GramianInverse',num2str(lausanneScaleBOLD),'Subject',num2str(N),'.mat']),'WcI_subj');
	A = SCvolnorm{N};
	A(~~eye(length(A))) = 0; % one of these matrices randomly had a NaN on diagonal
	A = NORMALIZE(A);
	subjectTransitionEnergy(N,:) = MIN_CONTROL_ENERGY(A,WcI_subj,Xo,Xf,T,false);
end

save(fullfile(savedir,['SubjectTransitionEnergy_k',num2str(numClusters),'.mat']),'subjectTransitionEnergy');