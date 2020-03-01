addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));

if scan == 'R'
    scanlab = {'Rest'};
elseif scan == 'N'
    scanlab = {'nBack'}; scanInd = scanInd * 0;
elseif scan == 'C'
    scanlab = {'RestComb','nBackComb'};

end
savedir = fullfile(masterdir,'/analyses/centroids');
mkdir(savedir);

for i = 1:numel(scanlab)
    
    load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat']);
    kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;
    clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
    clusterNamesUp = clusterAssignments.(['k',num2str(numClusters)]).clusterNamesUp;
    clusterNamesDown = clusterAssignments.(['k',num2str(numClusters)]).clusterNamesDown;
    kClusterCentroids = zeros(nparc,numClusters);
    
    for K = 1:numClusters
        kClusterCentroids(:,K) = mean(concTS(and(kClusterAssignments == K,scanInd == (i-1)),:),1)';
    end
    
    cd(savedir);

    save([scanlab{i},'ClusterCentroids_k',num2str(numClusters),name_root,'.mat'],'kClusterCentroids','clusterNames','clusterNamesDown','clusterNamesUp');

end

%compare rest vs. n-back
cd(savedir);

load([scanlab{1},'ClusterCentroids_k',num2str(numClusters),name_root,'.mat'],'kClusterCentroids');
restCentroids = kClusterCentroids;
load([scanlab{2},'ClusterCentroids_k',num2str(numClusters),name_root,'.mat'],'kClusterCentroids');
nBackCentroids = kClusterCentroids;

r = zeros(numClusters,1);
for K = 1:numClusters
    r(K) = corr(restCentroids(:,K),nBackCentroids(:,K));
end

% Fig S4b
%barcolors = [53 183 121; 68 1 84] / 255; 
RNcolors = {'458AC6','FFB75E'}; 
barcolors = hex2rgb(RNcolors);
f = figure; bar(r,'FaceColor',barcolors(1,:),'FaceAlpha',0.5); prettifyEJC
xticklabels(clusterNames); xtickangle(90); title({'Rest vs. n-back',['\mu_{r} = ',num2str(round(mean(r),2,'significant')),', \sigma_{r} = ',num2str(round(std(r),2,'significant'))]}); 
ylabel('Correlation');
COLOR_TICK_LABELS(true,false,numClusters);
f.PaperUnits = 'inches';
f.PaperSize = [3 1.5];
f.PaperPosition = [0 0 3 1.5];
saveas(f,['RvNCentroidSpatialCorr_k',num2str(numClusters),name_root,'.pdf'],'pdf');

% Get overall centroids

kClusterCentroids = zeros(nparc,numClusters);

for K = 1:numClusters
    kClusterCentroids(:,K) = mean(concTS(kClusterAssignments == K,:),1)';
end

[clusterNames,clusterReorder,clusterNamesSort] = NAME_CLUSTERS_ANGLE(kClusterCentroids);

clusterCorr = corr(kClusterCentroids);
% Fig S4a
imagesc(clusterCorr); colormap('plasma');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
clusterColors = {'AB484F','591A23', 'AA709F','527183','7E874B'};
if numClusters == 5
    COLOR_TICK_LABELS(true,true,numClusters,clusterColors);
end
h=colorbar; caxis([-1 1]); h.Ticks = [-1 0 1]; h.TickLabels = [-1 0 1];
title('Spatial Correlation');
set(gca,'FontSize',8);
f.PaperUnits = 'inches';
f.PaperSize = [2.7 2.7];
f.PaperPosition = [0 0 2.7 2.7];

saveas(f,['OverallCentroidSpatialCorr_k',num2str(numClusters),name_root,'.pdf'],'pdf');

clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
clusterNamesUp = clusterAssignments.(['k',num2str(numClusters)]).clusterNamesUp;
clusterNamesDown = clusterAssignments.(['k',num2str(numClusters)]).clusterNamesDown;
save(['OverallClusterCentroids_k',num2str(numClusters),name_root,'.mat'],'kClusterCentroids','clusterNames','clusterNamesDown','clusterNamesUp');
save(['Fig2a__OveralClusterCentroids_k',num2str(numClusters),name_root,'.mat'],'kClusterCentroids','clusterNames','clusterNamesDown','clusterNamesUp');
% overall centroids spatial anticorrelation