addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(basedir,'data','Gramians');
mkdir(savedir);

% pre-compute inverse gramian for each subject in groups of subjects
load(fullfile(datadir,['VolNormSC',num2str(lausanneScaleBOLD),'.mat']));

n_per_split = ceil(nobs/nsplits);	% calculate subjects per job based on # of jobs (nsplits) and # of subjects (nobs)
subj_range = (1+n_per_split*(split-1)):(n_per_split*split);
subj_range(subj_range > nobs) = [];		% truncate subj_range to equal total number of subjects... only happens on last job

T = 1; %control horizon
for N = 
	tic
	disp(['Subject ',num2str(N)])
	A = SCvolnorm{N};
	Anorm = (A / max(eigs(A,1))) - eye(length(A));	% normalize A to max eig --> new max eig = 1 --> one dominant eigenmode stationary over time
	WcI_subj = GRAMIAN(Anorm,T,false);
	save(fullfile(savedir,['GramianInverse',num2str(lausanneScaleBOLD),'Subject',num2str(N),'.mat']),'WcI_subj');
	toc
end

