addpaths;
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','hcpLR');
mkdir(savedir);
nparc = 462;
rTR = round(405); nTR = 405;

load(fullfile(basedir,['data/HCP_XH_concTS_R',num2str(rTR),'N',num2str(nTR),'.mat']),'concTS');

%% cluster
nreps = 15;
[HCPpartition,HCPcentroids,sumd] = kmeans(concTS,numClusters,'Distance',distanceMethod,'Replicates',nreps);
HCPcentroids = HCPcentroids';

%% compare centroids to pnc
%
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
PNCcentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
PNCnames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
PNCvsHCP = corr(HCPcentroids,PNCcentroids);
[~,shuffleIdx] = max(PNCvsHCP,[],1);
HCPcentroidsPNCorder = HCPcentroids(:,shuffleIdx);
clusterNames = NAME_CLUSTERS_ANGLE(HCPcentroidsPNCorder);
[clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(HCPcentroidsPNCorder);
kClusterCentroids = HCPcentroidsPNCorder;	% use for general plotting script
save(fullfile(savedir,['HCP_XHcentroids_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']),...
    'HCPcentroidsPNCorder','clusterNames','clusterNamesUp','clusterNamesDown','distanceMethod','kClusterCentroids');

f = figure;
PNCvsHCPorder = PNCvsHCP(shuffleIdx,:);
imagesc(PNCvsHCPorder); colormap('plasma');
% label diagonal with r-values
r_labels = cellstr(num2str(round(diag(PNCvsHCPorder),2,'significant')))';
text((1:numClusters)-0.1,(1:numClusters),r_labels,'Color','k','FontSize',6)

ylabel('HCP'); xlabel('PNC'); axis square
yticks(1:numClusters); xticks(1:numClusters);
yticklabels(clusterNames); xticklabels(PNCnames);
xtickangle(90);
COLOR_TICK_LABELS(true,true,numClusters);
h=colorbar; caxis([-1 1]); h.Ticks = [-1 0 1]; h.TickLabels = [-1 0 1];
title('Spatial Correlation');
set(gca,'FontSize',8);
f.PaperUnits = 'inches';
f.PaperSize = [2.7 2.7];
f.PaperPosition = [0 0 2.7 2.7];

saveas(f,fullfile(savedir,['HCP_XHClusterSpatialCorr_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.pdf']),'pdf');
%% within HCP spatial corr

f = figure;
HCPSpatialCorr = corr(HCPcentroidsPNCorder);
imagesc(HCPSpatialCorr); colormap('plasma');
r_labels = cellstr(num2str(round(diag(corr(HCPcentroidsPNCorder)),2,'significant')))';
text((1:numClusters)-0.1,(1:numClusters),r_labels,'Color','k','FontSize',6)
ylabel('HCP'); xlabel('HCP'); axis square
yticks(1:numClusters); xticks(1:numClusters);
yticklabels(clusterNames); xticklabels(PNCnames);
xtickangle(90);
COLOR_TICK_LABELS(true,true,numClusters);
h=colorbar; caxis([-1 1]); h.Ticks = [-1 0 1]; h.TickLabels = [-1 0 1];
title('Spatial Correlation');
set(gca,'FontSize',8);
f.PaperUnits = 'inches';
f.PaperSize = [2.7 2.7];
f.PaperPosition = [0 0 2.7 2.7];

saveas(f,fullfile(savedir,['HCP_XHBetweenClusterSpatialCorr_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.pdf']),'pdf');

save(fullfile(savedir,['FigS6b-c__HCPCentroidSpatialCorr_k',num2str(numClusters),'.mat']),'HCPSpatialCorr','PNCvsHCPorder');

%% reorder centroids based on similarity to original centroids

[~,shuffleIdx] = sort(shuffleIdx);	% reorder to index indiv. cluster not overall

idx = false(length(HCPpartition),numClusters);
% rest cluster K needs to become comb cluster shuffleIdx(K)

HCPpartitionPNCorder = zeros(length(HCPpartition),1);
for K = 1:numClusters 
	idx(:,K) = HCPpartition == K;  % find each rest cluster
end

for K = 1:numClusters
	HCPpartitionPNCorder(idx(:,K)) = shuffleIdx(K);  % replace with cluster # most similar to comb clusters
end

nsubjs = size(concTS,1) / (rTR+nTR);		% number of subjects is total # of TRs divided by # of TRs per subjects
HCPsubjInd = [repelem(1:nsubjs,rTR),repelem(1:nsubjs,nTR)];
HCPscanInd = [false(nsubjs*rTR,1);true(nsubjs*nTR,1)];

save(fullfile(savedir,['HCP_XHpartition_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']),...
    'HCPpartitionPNCorder','nsubjs','HCPsubjInd','HCPscanInd');