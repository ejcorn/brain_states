addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(basedir,'data','Group_DLWNull_Gramians');
mkdir(savedir);

% pre-compute a bunch of null matrices based on group representative SC that preserve edge weight-length relationships
% then compute inverse gramian and save

n_per_split = ceil(nperms/nsplits);	% calculate subjects per job based on # of jobs (nsplits) and # of subjects (nobs)
perm_range = (1+n_per_split*(split-1)):(n_per_split*split);

% load group representative SC
load(fullfile(masterdir,'analyses','control_energy',['GroupRepresentativeSC_FA_Laus',num2str(lausanneScaleBOLD),'.mat']),'A','D');

for P = perm_range
	tic
	disp(['Permutation ',num2str(P)])
    A_DLWnull = DLW_NULL(A,D);
	save(fullfile(savedir,['DLWNull_FA_Laus',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'A_DLWnull');
	toc
end

