addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(basedir,'data','Group_DLWNull_Gramians');
mkdir(savedir);

% pre-compute a bunch of spatial null matrices based on group representative SC
% then compute inverse gramian and save

n_per_split = ceil(nperms/nsplits);	% calculate subjects per job based on # of jobs (nsplits) and # of subjects (nobs)
perm_range = (1+n_per_split*(split-1)):(n_per_split*split);

% load group representative SC
load(fullfile(masterdir,'analyses','control_energy',['GroupRepresentativeSC_Laus',num2str(lausanneScaleBOLD),'.mat']),'A','D');
Anorm = (A / max(eig(A))) - eye(length(A));	% normalize A to max eig --> new max eig = 1 --> one dominant eigenmode stationary over time

T = 1; %control horizon
nbins = 11; nrewire = 1e5; 	% parameters for spatial null model
for P = perm_range
	tic
	disp(['Permutation ',num2str(P)])
	% empirically, null is more conservative i.e. control energies closer to actual matrix if you normalize then nullify (why I'm using Anorm as input)
	[~,A_DLWnull] = fcn_preserve_degseq_lengthdist(Anorm,D,nbins,nrewire); 
    A_DLWnull = A_DLWnull + A_DLWnull' - eye(length(A));;	% ensure matrix is symmetric    
	WcI_DLWnull = GRAMIAN(A_DLWnull,T,false);	% false means don't normalize
	save(fullfile(savedir,['DLWNullGramianInverse',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'WcI_DLWnull','A_DLWnull');
	toc
end

