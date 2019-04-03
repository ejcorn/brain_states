a=clock;
rng(a(6));
addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(basedir,'data','BootstrapGroupSC');
mkdir(savedir);

% pre-compute a bunch of spatial null matrices based on group representative SC
% then compute inverse gramian and save

n_per_split = ceil(nperms/nsplits);	% calculate subjects per job based on # of jobs (nsplits) and # of subjects (nobs)
perm_range = (1+n_per_split*(split-1)):(n_per_split*split);

% load group representative SC

load(fullfile(datadir,['VolNormSC',num2str(lausanneScaleBOLD),'.mat']));
load(fullfile(datadir,['Lausanne',num2str(lausanneScaleBOLD),'DistanceMatrix.mat']),'D');

T = 1; %control horizon
for P = perm_range
	tic
	disp(['Permutation ',num2str(P)])
	rng(P);		% make sure you get a unique bootstrap every time
	bootsamp = randi(nobs,1,nobs);		% sample with replacement
	A = cat(3,SCvolnorm{bootsamp});		% concatenate bootstrapped sample to 3d matrix
	A(repmat(~~eye(nparc),[1 1 nobs])) = 0;	% ensure diagonal is 0
	A = A .* repmat(~eye(nparc),[1 1 nobs]);  % get rid of diagonal elements
	frac = 1;                       % frac = 1 means your average matrix with have same # of elements as mean subject
	rightnodes = 230;
	hemi = [ones(rightnodes,1);2*ones(nparc - rightnodes,1)];
	G_boot = fcn_distance_dependent_threshold(A,D,hemi,frac);    % generate binary mask

	% you'll need to decide how to weight G -- in the past I've done the following:
	A(A == 0) = nan;    % set A = 0 to nan
	W = nanmean(A,3);   % nanmean (take into account only nonzero values when computing mean edge weight)
	A = G_boot.*W;          % set G = 1 to corresponding mean edge weight
	A(isnan(A)) = 0;

	Anorm = NORMALIZE(A);	% normalize A to max eig --> new max eig = 1 --> one dominant eigenmode stationary over time
	WcI_boot = GRAMIAN(Anorm,T,false);	% false means don't normalize
	A_boot = A;
	% save Gramian and binary mask
	save(fullfile(savedir,['BootstrapGroupRepresentativeSC',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'WcI_boot','G_boot','bootsamp','A_boot');
	toc
end

