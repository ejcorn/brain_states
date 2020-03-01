addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy');
mkdir(savedir);

%% get coordinates of spatial embedding

load(fullfile(datadir,['Lausanne',num2str(lausanneScaleBOLD),'DistanceMatrix.mat']),'D');

%% generate representative group SC matrix

load(fullfile(datadir,['FA',num2str(lausanneScaleBOLD),'.mat']));

A = zeros(nparc,nparc,nobs);
for N = 1:nobs
	A(:,:,N) = FA{N};
end
clear FA
         % size of connectivity matrices
A = A .* repmat(~eye(nparc),[1 1 nobs]);  % get rid of diagonal elements
frac = 1;                       % frac = 1 means your average matrix with have same # of elements as mean subject
rightnodes = 230;
hemi = [ones(rightnodes,1);2*ones(nparc - rightnodes,1)];
G = fcn_distance_dependent_threshold(A,D,hemi,frac);    % generate binary mask
%nbins = ceil(sqrt(nparc*(nparc-1)/2));	% set nbins equal to number of 
%G_nparc = fcn_group_bins(A,D,hemi,nbins);    % generate binary mask

% you'll need to decide how to weight G -- in the past I've done the following:
A(A == 0) = nan;    % set A = 0 to nan
W = nanmean(A,3);   % nanmean (take into account only nonzero values when computing mean edge weight)
A = G.*W;          % set G = 1 to corresponding mean edge weight
A(isnan(A)) = 0;

save(fullfile(savedir,['GroupRepresentativeSC_FA_Laus',num2str(lausanneScaleBOLD),'.mat']),'A','D','G');