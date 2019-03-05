addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(basedir,'data');

% pre-compute inverse gramian for each subject
load(fullfile(datadir,['VolNormSC',num2str(lausanneScaleBOLD),'.mat']));

T = 1; %control horizon
WcI = zeros(nparc,nparc,nobs);
for N = 1:nobs
	disp(['Subject ',num2str(N)])
	A = SCvolnorm{N};
	Anorm = (A / max(eigs(A,1))) - eye(length(A));	% normalize A to max eig --> new max eig = 1 --> one dominant eigenmode stationary over time
	WcI(:,:,N) = GRAMIAN(Anorm,T,false);
end

save(fullfile(savedir,['GramianInverse',num2str(lausanneScaleBOLD),'.mat']));