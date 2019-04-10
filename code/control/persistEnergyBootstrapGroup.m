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

bootstrappedPersistenceEnergy = zeros(nperms,numClusters);
bootsamps = zeros(nperms,nobs);
for P = 1:nperms
	disp(['Perm ',num2str(P)])
	load(fullfile(datadir,'BootstrapGroupSC',['BootstrapGroupRepresentativeSC',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'WcI_boot','G_boot','A_boot','bootsamp');
	% A = cat(3,SCvolnorm{bootsamp});	% construct bootstrapped stack of SC matrices
	% A(A == 0) = nan;    % set A = 0 to nan
	% W = nanmean(A,3);   % nanmean (take into account only nonzero values when computing mean edge weight)
	% A_boot = G_boot .* W; % set G = 1 to corresponding mean edge weight
	% A_boot(isnan(A_boot)) = 0;
	bootsamps(P,:) = bootsamp;	% save this to bootstrap persistence probabilities/dwell times as well
	A_boot = NORMALIZE(A_boot);
	bootstrappedPersistenceEnergy(P,:) = MIN_CONTROL_ENERGY(A_boot,WcI_boot,Xo,Xf,T,false);	%normalize A_boot
end

save(fullfile(savedir,['BootstrapPersistenceEnergy_k',num2str(numClusters),'.mat']),'bootstrappedPersistenceEnergy','bootsamps');

%{

G_cell = cell(nperms,1);
WcI_cell = cell(nperms,1);
A_cell = cell(nperms,1);
boot_cell = cell(nperms,1);
for P = 1:nperms
	disp(['Perm ',num2str(P)])
	load(fullfile(datadir,'BootstrapGroupSC',['BootstrapGroupRepresentativeSC',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'WcI_boot','G_boot','A_boot','bootsamp');
	WcI_cell{P} = WcI_boot;
	A_cell{P} = A_boot;
	G_cell{P} = G_boot;
	boot_cell{P} = bootsamp;
end

boot_u = uniquecell(boot_cell);
G_u = uniquecell(G_cell);
WcI_i = uniquecell(WcI_cell);

boots = zeros(nperms,nobs);
for i = 1:nperms
	rng(i);
	boots(i,:) = randi(nobs,1,nobs);
end

%}