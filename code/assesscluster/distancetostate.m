% goal of this script is to show distance to state centroids relative to 
% state time series
addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));
load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
partition = clusterAssignments.(['k',num2str(numClusters)]).partition;
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

%% load BOLD data from sample

savedir = fullfile(masterdir,'analyses','choosing_k'); mkdir(savedir);

%% time series distance vs. clusters

nreps = 5;
rTR = 120; nTR = 225;
subj = randperm(nobs,1);

partition_data = partition(subjInd == subj);
ts_data = concTS(subjInd==subj,:);  % get time series from one subject
d_to_centroids = zeros(nTR+rTR,numClusters);   % get distance of each time point to centroids
for K =1:numClusters
    d_to_centroids(:,K) = corr(ts_data',kClusterCentroids(:,K));
end

f=figure;
for i = 1:numClusters
    a=subplot(numClusters,1,i);    
    plot(partition_data == i,'b');
    hold on;
    plot(1-d_to_centroids_data(:,i),'r');
    line([rTR rTR],get(a,'YLim'),'Color',[0 0 0],'LineStyle',':');
    title([clusterNames{i},' State']);
    ylabel('\bf{\it{r}}','Interpreter','tex');
    ylim([-1 1]);
    if i == numClusters
        xlabel('Time');
    end
end

f.Units = 'centimeters';
f.PaperSize = [18 18];
f.PaperPosition = [0 0 18 18];
set(0,'CurrentFigure',f);
saveas(f,fullfile(savedir,['DistanceToStateVsState_',num2str(numClusters),'.pdf']));
%print(fullfile(savedir,'DistanceToStateVsState'),'-dpng','-r400');
