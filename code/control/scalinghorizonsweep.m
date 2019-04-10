addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy','sweep');
mkdir(savedir);

n_per_split = ceil(nperms/nsplits);	% calculate subjects per job based on # of jobs (nsplits) and # of subjects (nobs)
perm_range = (1+n_per_split*(split-1)):(n_per_split*split);

%% load states

% empirical states
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;

%% normalize states

StateMagnitude = sqrt(sum(kClusterCentroids.^2,1));
Xo = kClusterCentroids ./ StateMagnitude;		% normalize all states
Xf = Xo;	% start and end at same place, i.e. persistence energy

%% for each bootstrapped sample, compute persistence energy for each state for a range of control horizons and normalizations

T_rng = [0 10.^[-1 0 1 2]]; %control horizon range
T_rng = linspace(1,10,10);
n_horizons = length(T_rng);
% range of normalization factor c, where A is being normalized to 1/(c + max(eig(A))) - eye(length(A))
% when c = 0, max eig will be 0 and system goes to eigenvector associated with max eig at long time horizons
% in order to be more numerically precise, we use the function GRAMIAN here (precomputed for T = 1 in precomputeBootstrapSC.m)
% when c > 0, max eig will be < 0 and system goes to 0 at long time horizons
% we can compute the gramian for matrices normalized in this way very quickly using GRAMIAN_FAST (only doable b/c we are controlling from all nodes)
c_rng = 0:10; 
n_normalizations = length(c_rng);

sweepPersistenceEnergy = zeros(n_horizons,n_normalizations,numClusters,n_per_split);
for P = 1:n_per_split	% only store arrays at size of number of perms you process at once, then later reload in order
	disp(['Perm ',num2str(perm_range(P))])
	load(fullfile(datadir,'BootstrapGroupSC',['BootstrapGroupRepresentativeSC',num2str(lausanneScaleBOLD),'Perm',num2str(perm_range(P)),'.mat']),'WcI_boot','A_boot','bootsamp');
	tic
	for T_i = 1:n_horizons
		T = T_rng(T_i);
		for c_i = 1:n_normalizations
			c = c_rng(c_i);
			A = NORMALIZE(A_boot,c);
			WcI = GRAMIAN_FAST(A,T);
			sweepPersistenceEnergy(T_i,c_i,:,P) = MIN_CONTROL_ENERGY(A,WcI,Xo,Xf,T,false);
		end
	end
	toc
end

save(fullfile(savedir,['SweepHorizonNormalizationPersistenceEnergy_k',num2str(numClusters),'Split',num2str(split),'.mat']),'sweepPersistenceEnergy');
save(fullfile(savedir,['SweepHorizonNormalizationValueRanges.mat']),'T_rng','c_rng');