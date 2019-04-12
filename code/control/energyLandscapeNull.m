% purpose of this script is to compare transitions between random states to transitions from random states
basedir=pwd; name_root = 'ScanCLaus250Z0final'; numClusters = 5; nsplits = 25;
%addpaths;
addpath(genpath(pwd)); 
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
basedir=pwd;
masterdir = fullfile(basedir,'results',name_root);
%savedir = fullfile(masterdir,'analyses','control_energy');
savedir = fullfile(masterdir,'control_energy');
mkdir(fullfile(savedir,'energylandscape'));
%% load time series

%concTS = csvread(['data/ConcTSCSV_',name_root,'.csv']);
%concTS = randn(303255,nparc);
%% load states

load(fullfile(masterdir,'transitionprobabilities',['OverallClusterCentroids_k',num2str(numClusters),name_root,'.mat']));
Xo = kClusterCentroids;

%% load structure, normalize, compute gramian
load(fullfile(savedir,['GroupRepresentativeSC_Laus',num2str(lausanneScaleBOLD),'.mat']),'A','D');
%load(fullfile(datadir,'BootstrapGroupSC',['BootstrapGroupRepresentativeSC',num2str(lausanneScaleBOLD),'Perm1.mat']),'WcI_boot','G_boot','A_boot','bootsamp');
c = 0;
A = NORMALIZE(A,c);

T = 1; %control horizon
WcI = GRAMIAN_FAST(A,T);	% compute gramian
%% select random states

n_states = 1000; nTRs = size(concTS,1);
shuffidx1 = randi(nTRs,1,n_states*1.5);
shuffidx2 = randi(nTRs,1,n_states*1.5);
% make sure you only have transitions, vs. staying in same states
transidx = shuffidx1 ~= shuffidx2;
shuffidx1 = shuffidx1(transidx); shuffidx1 = shuffidx1(1:n_states);
shuffidx2 = shuffidx2(transidx); shuffidx2 = shuffidx2(1:n_states);

%% compute energy from random states to random states -- centroids should be easy to reach

% select two sets of random BOLD time points as nulls
Xo_Null = concTS(shuffidx1,:)';
Xf_Null = concTS(shuffidx2,:)';

% normalize states
InterStateDistance = sqrt(sum((Xf_Null - Xo_Null).^2,1));
Xo_Null = Xo_Null./InterStateDistance;
Xf_Null = Xf_Null./InterStateDistance;
disp('After normalizing states, distance between every pair should be 1:');
disp(sum((Xo_Null - Xf_Null).^2,1));
Random_to_random = MIN_CONTROL_ENERGY(A,WcI,Xo_Null,Xf_Null,T,false);

Random_to_centroids = zeros(n_states,numClusters);

for K = 1:numClusters
    % repeat each centroid n_states times to use as final state
    Xf = repmat(kClusterCentroids(:,K),1,n_states);
    Xo_Null = concTS(shuffidx1,:)';
    InterStateDistance = sqrt(sum((Xf-Xo_Null).^2,1));
    Xo_Null = Xo_Null./InterStateDistance;
    Xf = Xf./InterStateDistance;
    Random_to_centroids(:,K) = MIN_CONTROL_ENERGY(A,WcI,Xo_Null,Xf,T,false);
end

%% plot
f=figure; 
histogram(Random_to_random'); 
hold on;
pvals = zeros(numClusters,1);
labels = ['BOLD';clusterNames];
for K = 1:numClusters
    histogram(Random_to_centroids(:,K));
    % Bonferroni-corrected one-tailed non-parametric p-value for energy of random to centroids >
    % random to random
    pvals(K) = numClusters*mean(Random_to_centroids(:,K) > Random_to_random');
    if pvals(K) == 0
        pval_label = ['p < ',num2str(1/n_states)];
    else
        pval_label = ['p = ',num2str(round(pvals(K),2,'significant'))];
    end
    labels{K+1} = [labels{K+1},', ',pval_label];
end
xlabel('Min. Control Energy');
ylabel('Count');
title({'BOLD Frames to BOLD Frames vs.',' BOLD Frames to Centroids'});
legend(labels)
set(gca,'FontSize',8);
f.PaperUnits = 'centimeters';
f.PaperSize = [9 9];
f.PaperPosition = [0 0 9 9];
saveas(f,fullfile(savedir,'energylandscape',['MinControlEnergyBOLDtoBOLDVsBOLDtoCentroids_k',num2str(numClusters),'c',num2str(c),...
    'T',num2str(T),'.pdf']));

%% centroids to random -- centroids should be hard to leave

Centroids_to_random = zeros(n_states,numClusters);

for K = 1:numClusters
    % repeat each centroid n_states times to use as initial state
    Xo = repmat(kClusterCentroids(:,K),1,n_states);
    % use random BOLD time points as final state
    Xf_Null = concTS(shuffidx1,:)';
    InterStateDistance = sqrt(sum((Xf_Null - Xo).^2,1));
    Xf_Null = Xf_Null./InterStateDistance;
    Xo = Xo./InterStateDistance;
    Centroids_to_random(:,K) = MIN_CONTROL_ENERGY(A,WcI,Xo,Xf_Null,T,false);
end

%% plot
f=figure; 
histogram(Random_to_random'); 
hold on;
pvals = zeros(numClusters,1);
labels = ['BOLD';clusterNames];
for K = 1:numClusters
    histogram(Centroids_to_random(:,K));
    % one-tailed non-parametric p-value for energy of centroids to random <
    % random to random
    pvals(K) = numClusters*mean(Centroids_to_random(:,K) < Random_to_random');
    if pvals(K) == 0
        pval_label = ['p < ',num2str(1/n_states)];
    else
        pval_label = ['p = ',num2str(round(pvals(K),2,'significant'))];
    end
    labels{K+1} = [labels{K+1},', ',pval_label];
end
xlabel('Min. Control Energy');
ylabel('Count');
title({'BOLD Frames to BOLD Frames vs.','Centroids to BOLD Frames'});
legend(labels,'Location','northwest')
set(gca,'FontSize',8);
f.PaperUnits = 'centimeters';
f.PaperSize = [9 9];
f.PaperPosition = [0 0 9 9];
saveas(f,fullfile(savedir,'energylandscape',['MinControlEnergyBOLDtoBOLDVsCentroidstoBold_k',num2str(numClusters),'c',num2str(c),...
    'T',num2str(T),'.pdf']));