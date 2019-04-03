% compare uniform phase-randomized null vs. actual data for rest and n-back
% separately then together
addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
savedir = fullfile(masterdir,'analyses','choosing_k');
mkdir(savedir);

%% load data 
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));
load(fullfile(datadir,['iprtimeseries',name_root,'.mat']),'iprTS');
load(fullfile(datadir,['randtimeseries',name_root,'.mat']),'randTS');

load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
partition_data = clusterAssignments.(['k',num2str(numClusters)]).partition;

load(fullfile(masterdir,'clusterAssignments',['randnpartition_k',num2str(numClusters)]),'partition_randn')
load(fullfile(masterdir,'clusterAssignments',['IPRpartition_k',num2str(numClusters)]),'partition_ipr')

%%
n_sil = 40;
sil_mask = ismember(subjInd,randperm(nobs,n_sil));

f = figure; 
subplot(1,4,1);
[S_randn,~] = silhouette(randTS(sil_mask,:),partition_randn(sil_mask),distanceMethod);
title({'Independent','Random Gaussians'});
set(gca,'FontSize',8);

subplot(1,4,2);
[S_ipr,~] = silhouette(iprTS(sil_mask,:),partition_ipr(sil_mask),distanceMethod);
title({'Autocorrelation-Preserving','Null'});
set(gca,'FontSize',8);

subplot(1,4,3);
[S_data,~] = silhouette(concTS(sil_mask,:),partition_data(sil_mask),distanceMethod);
title('Real Data');
set(gca,'FontSize',8);

subplot(1,4,4);
histogram(S_data,'EdgeAlpha',0.3); hold on;
histogram(S_ipr,'EdgeAlpha',0.3);
xlim([min([S_data;S_ipr]) max([S_data;S_ipr])]);
xlabel('Silhouette Value');
ylabel('# of TRs');
[~,p,ci] = ttest2(S_data,S_ipr);
title(['\mu_{real}-\mu_{ipr} = ',num2str(round(mean(ci),2,'significant')),', p = ',num2str(round(p,2,'significant'))])
set(gca,'FontSize',8);

f.PaperUnits = 'centimeters';
f.PaperSize = [21 4];
f.PaperPosition = [0 0 21 4];
set(0,'CurrentFigure',f);
saveas(f,fullfile(savedir,['SilhouetteVsNulls_k',num2str(numClusters),'.pdf']));
%print(fullfile(savedir,'SilhouetteVsNulls.png'),'-dpng','-r400');

%%

f = figure;
histogram(S_ipr); hold on;
histogram(S_data);
xlabel('Silhouette Value');
ylabel('Count');
[~,p,ci,stats] = ttest2(S_data,S_ipr);
stats.tstat
stats.df
title(['\mu_{real}-\mu_{ipr} = ',num2str(round(mean(ci),2,'significant')),', p = ',num2str(round(p,2,'significant'))])
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 8 8];
f.PaperSize = [8 8];
saveas(f,fullfile(savedir,['SilhouetteRealvsIPRHistogram.pdf']));