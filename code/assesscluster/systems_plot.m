addpaths;
masterdir = fullfile(basedir,'results',name_root);

%%
load(fullfile(masterdir,'clusterAssignments',['k',num2str(numClusters),name_root,'.mat']));
overallPartition = clusterAssignments.(['k',num2str(numClusters)]).partition;
centroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
overallNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
savedir = fullfile(masterdir,'analyses','centroids');
mkdir(savedir);

%%
[nparc,numClusters] = size(centroids);

[~,~,~,net7angle] = NAME_CLUSTERS_ANGLE(centroids);

YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-'};

%% plot

clusterColors = GET_CLUSTER_COLORS(numClusters);
%YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;247 253 205; 218 166 86; 199 109 117] / 255;
YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;250 218 94; 218 166 86; 199 109 117] / 255;
YeoColors = [YeoColors;YeoColors];

f = figure;
imagesc(net7angle); ax = gca;
colormap('plasma')
set(ax,'xaxisLocation','top')
xticks(1:14); xticklabels(YeoNetNames);
xtickangle(90);
for K = 1:14
	ax.XTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
	YeoColors(K,:), ax.XTickLabel{K});
end
yticks(1:numClusters); COLOR_TICK_LABELS(false,true,numClusters);
set(ax,'FontSize',8);

%%

[~,~,net7angle_Up,net7angle_Down] = NAME_CLUSTERS_UP_DOWN(centroids);
YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};
numNets = numel(YeoNetNames);

%% plot

sim = [net7angle_Up;net7angle_Down];
reOrder = [1 6 2 7 3 8 4 9 5 10];   % order 1+ 1- 2+ 2- etc.
clusterColors = {'AB484F','591A23', 'AA709F','527183','7E874B'};
clusterColors(6:10) = clusterColors;
clusterColors = hex2rgb(clusterColors);
sim = sim(reOrder,:);
clusterColors = clusterColors(reOrder,:);
yticklabs = {'1+','2+','3+','4+','5+','1-','2-','3-','4-','5-'};
f = figure;
imagesc(sim); ax = gca;
colormap('plasma')
set(ax,'xaxisLocation','top')
xticks(1:(numNets*2)); xticklabels(YeoNetNames);
xtickangle(90);
for K = 1:(numNets*2)
	ax.XTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
	YeoColors(K,:), ax.XTickLabel{K});
end
yticks(1:numClusters*2); yticklabels(yticklabs(reOrder)');
for K = 1:(numClusters*2)
	ax.YTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
	clusterColors(K,:), ax.YTickLabel{K});
end
set(ax,'FontSize',8);
set(ax,'TickLength',[0 0]);
box off

for i = 2:2:8
    line([0 7.5],[i + 0.5 i + 0.5],'LineWidth',1,...
        'color',[0.5 0.5 0.5]);
end
[~,sysIdx] = max(sim,[],2);
for i = 1:length(sysIdx)
    rectangle('Position',[sysIdx(i)-0.5,i-0.5,1,1],'EdgeColor','k',...
        'LineWidth',2);
end

f.PaperUnits = 'inches';
f.PaperSize = [2.5 2];
f.PaperPosition = [0 0 2.5 2];
saveas(f,fullfile(savedir,['Systems_k',num2str(numClusters),name_root,'.pdf']));


%% make radial plots

clusterColors = GET_CLUSTER_COLORS(numClusters);

clusterColors = hex2rgb(clusterColors);
netAngle = linspace(0,2*pi,numNets+1);
thetaNames = YeoNetNames; thetaNames{8} = '';
f=figure;
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[net7angle_Up(K,:) net7angle_Up(K,1)],'k');
    polarplot(netAngle,[net7angle_Down(K,:) net7angle_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
        YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',6);
    title(overallNames{K},'Color',clusterColors(K,:),'FontSize',8);
end
f.PaperUnits = 'inches';
f.PaperSize = [8 1.5];
f.PaperPosition = [0 0 8 1.5];
saveas(f,fullfile(savedir,['SystemsRadial_k',num2str(numClusters),name_root,'.pdf']));

% for source data file
save(fullfile(savedir,['Fig2b__YeoSystemAlignment_k',num2str(numClusters),'.mat']),'netAngle','net7angle_Up','net7angle_Down');
