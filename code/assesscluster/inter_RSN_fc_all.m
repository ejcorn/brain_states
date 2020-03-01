% Makes Fig S8a

addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));
load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
partition_data = clusterAssignments.(['k',num2str(numClusters)]).partition;
load(['data/yeo7netlabelsLaus',num2str(lausanneScaleBOLD),'.mat'])

savedir = fullfile(masterdir,'analyses','fc'); mkdir(savedir);

network7labels = network7labels(1:(end-1)); % remove brainstem
[netlabels_sort,I] = sort(network7labels);
concTS = concTS(:,I);

%% select rest or n-back
names = {'rest','n-back'};
r_n = [0 1]; %0 for rest, 1 for n-back
r_n_name = [names{r_n+1}];

%% get functional connectivity

fc = zeros(nparc,nparc,nobs,length(r_n));	% separately for rest and n-back, if using both

for scan = 1:length(r_n)
	for i = 1:nobs
	    fc_i = corr(concTS(subjInd' == i & ismember(scanInd,r_n(scan)),:));
	    fc_i(logical(eye(nparc))) = NaN;	% make diagonal nan to exclude with nanmean
	    fc(:,:,i,scan) = fc_i;
	end
end

fc = mean(fc,4);	% average each connection for each subject across rest and n-back if using both

%% compare FC within & between systems
% goal
% DMN- state: DMN down, VIS up, SOM up
% is this pattern expected from functional connectivity?
% cofluctuation of regions within VIS, SOM, DMN separately is expected from
% temporal correlation structure, i.e. because of within network connectivity
% but why this combo of cofluctuating networks?
% that reflects between network connectivity or off-diagonal structure

%%

numNets = 7;
SystemAverageConnectivity = zeros(numNets);
SystemAverageConnectivity_pval = zeros(numNets);
SystemAverageConnectivity_tstat = zeros(numNets);
r_labels = cell(numNets);
for Net1 = 1:numNets
	Net1mask = netlabels_sort == Net1;
	for Net2 = 1:numNets
		Net2mask = netlabels_sort == Net2;
		NetPair_FC = fc(Net1mask,Net2mask,:);		% extract connections between Network1 and Network2
		meanNetPair_FC = squeeze(nanmean(nanmean(NetPair_FC,1),2)); % average correlations
		SystemAverageConnectivity(Net1,Net2) = mean(meanNetPair_FC);
		r_labels{Net1,Net2} = num2str(round(mean(meanNetPair_FC),2));
		[~,SystemAverageConnectivity_pval(Net1,Net2),~,stats] = ttest(meanNetPair_FC);
		SystemAverageConnectivity_tstat(Net1,Net2) = stats.tstat;
	end
end

[yp,xp] = find((SystemAverageConnectivity_pval < 0.05));
[yr,xr] = find((ones(size(SystemAverageConnectivity))));

netNames = {'VIS','SOM','DAN','VAN','LIM','FPN','DMN'};

SystemAverageConnectivity	% show values to write in
SystemAverageConnectivity_tstat
SystemAverageConnectivity_pval

% Fig S8a
f=figure;
imagesc(SystemAverageConnectivity); 
text(xp,yp,'*','Color','w');
text(xr-0.25,yr+0.2,reshape(r_labels,1,[]),'Color','k','FontSize',6)
xticks(1:numNets); xticklabels(netNames); xtickangle(90);
yticks(1:numNets); yticklabels(netNames);
maxval = max(max(abs(SystemAverageConnectivity)));
colorbar; caxis([-maxval maxval]); colormap('plasma');

f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 9 6];
f.PaperSize = [9 6];
saveas(f,fullfile(savedir,[r_n_name,'MeanInterYeoSystemFC.pdf']));
