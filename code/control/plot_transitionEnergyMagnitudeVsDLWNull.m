addpaths;
%%
% to run locally:
%{
% switch all basedir to basedir local
basedir_local = '~/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states';
addpath(genpath(fullfile(basedir_local,'code')));
addpath(genpath('~/Dropbox/Cornblath_Bassett_Projects/code/brainmapping2'));
c=0; T = 2; numClusters = 5;
name_root = 'ScanCLaus250Z0final';
%}
%%
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
%masterdir = fullfile(basedir_local,'results',name_root,'control_energy');
masterdir = fullfile(basedir,'results',name_root,'analyses','control_energy'); % for CFN
savedir = fullfile(masterdir,'null_network_energies');
mkdir(savedir);

% load transition energies for null
load(fullfile(masterdir,['DLWNullsTransitionEnergy_FA_k',num2str(numClusters),'c',num2str(c),'T',num2str(T),'.mat']),'DLWNullTransitionEnergy');

% load transition energies for group representative matrix
load(fullfile(masterdir,'Tsweep',['c',num2str(c),'_k',num2str(numClusters)],...
    ['GroupAverageTransitionEnergies_k',num2str(numClusters),'.mat']));

% load cluster names
load(['results/',name_root,'/analyses/centroids/OverallClusterCentroids_k',num2str(numClusters),name_root,'.mat'],'clusterNames');

%% compute pvalues
[~,T_idx] = min(abs(T_rng-T)); % find index of T closest to desired T (will be the T_min identified at the end of transitionEnergyDynamicsGroupv2_TSweep.m)
E_matrix = reshape(E_full(:,T_idx)',[numClusters numClusters])';

% reshape BCT null into matrix
[num_transitions,nperms,~] = size(E_BCTnull);
BCTNullTransitionEnergyMatrix = zeros(numClusters,numClusters,nperms);

for P = 1:nperms
	BCTNullTransitionEnergyMatrix(:,:,P) = reshape(E_BCTnull(:,P,T_idx),[numClusters numClusters])';
end

% reshape DLW null into matrix
[nperms,num_transitions] = size(DLWNullTransitionEnergy);
DLWNullTransitionEnergyMatrix = zeros(numClusters,numClusters,nperms);

for P = 1:nperms
    DLWNullTransitionEnergyMatrix(:,:,P) = reshape(DLWNullTransitionEnergy(P,:),[numClusters numClusters])';
end

sig_thresh = 0.05 / num_transitions; % bonferroni correct over all transitions
pvals_onetail_DLW = mean(E_matrix > DLWNullTransitionEnergyMatrix,3); % get p-value as % of times real energy is higher than energy in null matrix for a transition
pvals_onetail_BCT = mean(E_matrix > BCTNullTransitionEnergyMatrix,3); % get p-value as % of times real energy is higher than energy in null matrix for a transition

%% plot energy matrix
f=figure;
imagesc(E_matrix); 
h = colorbar; 
ylabel(h,'E_{min}');
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
[y,x] = find(pvals_onetail_BCT.*~eye(numClusters) < sig_thresh);
text(x-.12,y,'*','Color',[0.5 0.5 0.5]);
[y,x] = find(pvals_onetail_DLW.*~eye(numClusters) < sig_thresh);
text(x-.12,y+.28,'*','Color','w');
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Initial State'); xlabel('Final State');
title('Minimum Control Energy');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 6 6];
f.PaperSize = [6 6];

save(fullfile(savedir,'Fig5b__ControlEnergyMagnitudeVsNulls.mat'),'pvals_onetail_DLW','pvals_onetail_BCT','E_matrix');
saveas(f,fullfile(savedir,['TransitionEnergyVsDLWNull_c',num2str(c),'T',num2str(T),'_k',num2str(numClusters),'.pdf']));