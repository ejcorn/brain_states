clear all; close all;clc
basedir = '~/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states/';
cd(basedir);
addpath(genpath('code'))
%% set inputs

savedir = fullfile(basedir,'results','example');mkdir(savedir);		% set save directory
distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
numClusters = 5; % set number of clusters -- this requires a separate process of your choosing to select the best number of clusters. See https://arxiv.org/abs/1007.1075 for a useful discussion. 
nreps = 20;	% how many times to repeat clustering. will choose lowest error solution

%% load BOLD data

% replace TS with your BOLD data formatted as a T-by-nparc matrix
% where T is the number of time points and nparc is the number of ROIs
D = load('data/HCP_XH_concTS_R405N405.mat'); % HCP U100 sample in lausanne250 parcellation
TS = D.concTS;
rTR = round(405); nTR = 405; % use sa
nsubjs = size(TS,1) / (rTR+nTR);		% number of subjects is total # of TRs divided by # of TRs per subjects
subjInd = [repelem(1:nsubjs,rTR),repelem(1:nsubjs,nTR)]'; % index data from each subject
scanInd = [false(nsubjs*rTR,1);true(nsubjs*nTR,1)]; % index scans within your time series data
scanlab = {'rest','n-back'}; % text labels for scan indices
[T,nparc] = size(TS);

%% identify brain states
[partition,~,sumd] = kmeans(TS,numClusters,'Distance', distanceMethod,'Replicates',nreps);

% see supplement for information on evaluating choices for the number of clusters
% in our experience, 4-6 clusters are around the elbow of the variance
% explained curve and are robust to subsampling.
%% compute centroids and plot
centroids = GET_CENTROIDS(TS,partition,numClusters);
% name clusters based on alignment with Yeo resting state networks
clusterNames = NAME_CLUSTERS_ANGLE(centroids);  % need to add a prior Yeo partition labels for your parcellation
[clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(centroids);  % need to add a prior Yeo partition labels for your parcellation

f = figure;
subplot(1,2,1); imagesc(centroids); title('Centroids'); xticks(1:numClusters); xticklabels(clusterNames);
colormap('plasma'); axis square; colorbar; set(gca,'FontSize',8); COLOR_TICK_LABELS(true,false,numClusters);
subplot(1,2,2); imagesc(corr(centroids)); title('Centroid Similarity'); colorbar; caxis([-1 1]); 
colormap('plasma'); axis square; set(gca,'FontSize',8); xticks(1:numClusters); yticks(1:numClusters); 
xticklabels(clusterNames); yticklabels(clusterNames); xtickangle(90);
COLOR_TICK_LABELS(true,true,numClusters);
f.PaperUnits = 'inches';
f.PaperSize = [4 2];
f.PaperPosition = [0 0 4 2];
saveas(f,fullfile(savedir,['Centroids_k',num2str(numClusters),'.pdf']));

% can use python scripts to make surface plots
%% transition probabilities

% transition probability matrices:
% restTransitionProbabilityMats: each element reflects probability that
% state j occurs at the TR after state i, given that you are in state i (and not in the last TR of a scan)
% restTransitionProbabilityMatsNoPersist: each element reflectts
% probability that state j follows state i in the next state change,
% given that you are in state i (and not in the last TR of a scan). This
% *NoPersist version is used in the manuscript but both are interesting.

% get resting state transition probability matrices -- 2D is flattened by row for
% regressions
[restTransitionProbability2D,restTransitionProbabilityMats] = GET_TRANS_PROBS(partition(scanInd == 0),subjInd(scanInd == 0));
[restTransitionProbabilityNoPersist2D,restTransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(partition(scanInd == 0),subjInd(scanInd == 0));

% get transition probabilities within two back blocks
load data/HCP_TwoBackBlock.mat % load indicator vector specifying location of two back blocks

twobackTransitionProbabilityNoPersist2D = zeros(nsubjs,numClusters^2);
twobackTransitionProbabilityMatsNoPersist = zeros(nsubjs,numClusters,numClusters);	% preallocate transition probability matrices
for N = 1:nsubjs
	subjPartition = partition(subjInd == N & scanInd == 1);
	[twobackTransitionProbabilityMatsNoPersist(N,:,:),twobackTransitionProbabilityNoPersist2D(N,:)] = GET_BLOCK_TRANS_PROBS_NO_PERSIST(subjPartition,TwoBackBlock,numClusters);
end
grpAvgRest = squeeze(mean(restTransitionProbabilityMatsNoPersist,1)) .* ~eye(numClusters);
% nans occur for transitions from states that are not present at all for a subject or within the tested blocks. 
% transitions to that state are 0
grpAvg2Back = squeeze(nanmean(twobackTransitionProbabilityMatsNoPersist,1)) .* ~eye(numClusters); 

%% permutation testing to compare transition probability matrices
disp('start permutation testing')
nperms = 100000;
pvals_twotail = PERM_TEST(twobackTransitionProbabilityMatsNoPersist,restTransitionProbabilityMatsNoPersist,nperms);

%% plot transition probabilities
maxVal = max(max([grpAvgRest,grpAvg2Back])); % sync color scales

f = figure;

subplot(1,3,1);
imagesc(grpAvgRest);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next New State');
title('Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,3,2);
imagesc(grpAvg2Back);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);

ylabel('Current State'); xlabel('Next New State');
title('2-back');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,3,3);
nBackMinusRestTP = (grpAvg2Back-grpAvgRest);
imagesc(nBackMinusRestTP.*~eye(numClusters)); colormap('plasma');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Current State'); xlabel('Next New State');
sig_thresh = 0.05 / numClusters^2;      % bonferroni correction, for two-tailed p-values so only
[y,x] = find(pvals_twotail.*~eye(numClusters) < sig_thresh);
text(x-.12,y+.12,'*','Color','w');
caxis_bound = max(max(abs(nBackMinusRestTP.*~eye(numClusters))));
h = colorbar; ylabel(h,'2-back - rest'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('n-back > Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
% plot transition probability matrices
saveas(f,fullfile(savedir,['Rest2BackTransProbs_k',num2str(numClusters),'.pdf']),'pdf');

% see code/transprobs/plot_transprob_digraph.m for visualizing TPs as a
% directed network

%% control energy

load data/ExampleSC_FA_Laus250.mat % load example group average structural A matrix -- this is from PNC while fMRI is HCP so DTI is younger than fMRI here
c = 0; T = 5; % set time scale parameters based on values from paper
Anorm = NORMALIZE(A,c); % normalize A by maximum eigenvalue - eye(N) to make marginally stable

% define x0 and xf, initial and final states as cluster centroids for each
% state transition
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

x0 = centroids(:,Xo_ind);
xf = centroids(:,Xf_ind); % now each column of x0 and xf represent state transitions
WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
E_full = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition

% compute weighted control energy:
load(fullfile(basedir,'data','yeo7netlabelsLaus250.mat'));
InputVector = ismember(network7labels(1:nparc),1); % weight input towards visual system as in task
B = InputVector .*eye(nparc) + eye(nparc); % construct input matrix allowing input only into selected regions in InputVector
E_weighted = zeros(1,numClusters^2);
for transition = 1:numClusters^2    
    [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
    E_weighted(transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
end
