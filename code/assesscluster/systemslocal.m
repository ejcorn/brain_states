addpaths;
masterdir = fullfile(basedir,'results',name_root);

%%
load(fullfile(masterdir,'clusterAssignments',['k',num2str(numClusters),name_root,'.mat']));
overallPartition = clusterAssignments.(['k',num2str(numClusters)]).partition;
centroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
overallNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
savedir = fullfile(masterdir,'assign');
mkdir(savedir);
cd(savedir);
%%
[nparc,numClusters] = size(centroids);

if nparc > 400
    load('data/yeo7netlabelsLaus250.mat'); network7labels = network7labels(1:nparc);
else
    load('data/yeo7netlabelsLaus125.mat'); network7labels = network7labels(1:nparc);
end

binaryNetVectors = zeros(nparc,14);

for I = 1:7
    binaryNetVectors(:,I) = double(network7labels == I);
    binaryNetVectors(:,I+7) = -double(network7labels == I);
end

YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-'};

% calculate E.D. from binary state vector to centroids

net7ED = zeros(numClusters,14);

for K = 1:numClusters
    for B = 1:14
        net7ED(K,B) = dot(centroids(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

%% plot

clusterColors = {'AB484F','591A23', 'AA709F','527183','7E874B'};
%YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;247 253 205; 218 166 86; 199 109 117] / 255;
YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;250 218 94; 218 166 86; 199 109 117] / 255;
YeoColors = [YeoColors;YeoColors];

f = figure;
imagesc(net7ED); ax = gca;
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

numNets = 7;
binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == network7labels);

YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};

% calculate E.D. from binary state vector to centroids

centroids_up = centroids .* (centroids > 0);
centroids_down = -1 * centroids .* (centroids < 0);     % make negative activity positive and get rid of positive activity

net7ED_Up = zeros(numClusters,numNets);
net7ED_Down = zeros(numClusters,numNets);

for K = 1:numClusters
    for B = 1:numNets
        net7ED_Up(K,B) = dot(centroids_up(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
        net7ED_Down(K,B) = dot(centroids_down(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

%% plot
sim = [net7ED_Up;net7ED_Down];
reOrder = [1 6 2 7 3 8 4 9 5 10];
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
xticks(1:14); xticklabels(YeoNetNames);
xtickangle(90);
for K = 1:14
	ax.XTickLabel{K} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
	YeoColors(K,:), ax.XTickLabel{K});
end
yticks(1:numClusters*2); yticklabels(yticklabs(reOrder)');
for K = 1:10
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
saveas(f,['Systems_k',num2str(numClusters),name_root,'.pdf'])

%% make radial plots

clusterColors = {'AB484F','591A23', 'AA709F','527183','7E874B'};
clusterColors = hex2rgb(clusterColors);
netAngle = linspace(0,2*pi,numNets+1);
thetaNames = YeoNetNames; thetaNames{8} = '';
f=figure;
for K = 1:numClusters
    ax = subplot(1,5,K,polaraxes); hold on
    polarplot(netAngle,[net7ED_Up(K,:) net7ED_Up(K,1)],'k');
    polarplot(netAngle,[net7ED_Down(K,:) net7ED_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rticks([0.3 0.6]); rticklabels({'0.3','0.6'});
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
saveas(f,['SystemsRadial_k',num2str(numClusters),name_root,'.pdf'])
