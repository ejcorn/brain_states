% makes FigS8

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
DMN = 7; VIS = 1; SOM = 2;
DMN_mask = netlabels_sort == DMN;
VIS_mask = netlabels_sort == VIS;
SOM_mask = netlabels_sort == SOM;

% compute correlations within a priori systems
triu_mask = repmat(~tril(true(sum(DMN_mask))),[1 1 nobs]);
within_DMN = fc(DMN_mask,DMN_mask,:);
within_DMN = squeeze(nanmean(nanmean(within_DMN,1),2)); % average corr
%within_DMN = within_DMN(triu_mask);  % all corrs

triu_mask = repmat(~tril(true(sum(VIS_mask))),[1 1 nobs]);
within_VIS = fc(VIS_mask,VIS_mask,:);
within_VIS = squeeze(nanmean(nanmean(within_VIS,1),2)); % average corr
%within_VIS = within_VIS(triu_mask);  % all corrs

triu_mask = repmat(~tril(true(sum(SOM_mask))),[1 1 nobs]);
within_SOM = fc(SOM_mask,SOM_mask,:);
within_SOM = squeeze(nanmean(nanmean(within_SOM,1),2)); % average corr
%within_SOM = within_SOM(triu_mask);  % all corrs

numbins = 10;

% Fig S8b-d
f=figure;
subplot(1,3,1);
histogram(within_SOM,numbins,'EdgeAlpha',0.5); hold on;
histogram(within_VIS,numbins,'EdgeAlpha',0.5);
histogram(within_DMN,numbins,'EdgeAlpha',0.5);
xlim([-0.8 0.8])
%legend({'DMN','VIS','SOM'})
title('Within System');
ylabel('# of Subjects');
xlabel('r');
%% between network connectivity

% compute correlations within a priori systems
between_DMN_VIS = fc(DMN_mask,VIS_mask,:);
between_DMN_VIS = squeeze(nanmean(nanmean(between_DMN_VIS,1),2)); % average corr

between_DMN_SOM = fc(DMN_mask,SOM_mask,:);
between_DMN_SOM = squeeze(nanmean(nanmean(between_DMN_SOM,1),2)); % average corr

between_VIS_SOM = fc(VIS_mask,SOM_mask,:);
between_VIS_SOM = squeeze(nanmean(nanmean(between_VIS_SOM,1),2)); % average corr

subplot(1,3,2);
histogram(between_VIS_SOM,numbins,'EdgeAlpha',0.5); hold on;
histogram(between_DMN_SOM,numbins,'EdgeAlpha',0.5);
histogram(between_DMN_VIS,numbins,'EdgeAlpha',0.5);
xlim([-0.8 0.8])
%legend({'DMN-VIS','DMN-SOM','VIS-SOM'},'Location','northwest')
title('Between System');
ylabel('# of Subjects');
xlabel('r');

%% 
% ***show the distribution of correlations between time points in DMN-
% with binary DMN, VIS, SOM

DMN_low = find(strcmp(clusterNames,'DMN-'));
Yeo_ref = double(netlabels_sort == SOM) + double(netlabels_sort == VIS) - double(netlabels_sort == DMN);
% correlate across only regions in DMN, VIS and SOM
DMN_minus = corr(concTS(partition_data == DMN_low,Yeo_ref ~= 0)',Yeo_ref(Yeo_ref ~= 0));

subplot(1,3,3);
histogram(DMN_minus,numbins,'EdgeAlpha',0.5);
ylabel('# of Frames');
xlabel('r');
xlim([-0.8 0.8]);
title({'TRs in DMN- with','DMN_{low}, SOM_{high}, VIS_{high}'});
%f.Renderer='Painters';
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 21 6];
f.PaperSize = [21 6];
saveas(f,fullfile(savedir,[r_n_name,'FCWithinBetween.pdf']));

%% export stats for latex table
r_tab = zeros(7,4);
[~,p,ci] = ttest(within_DMN);
r_tab(1,:) = [mean(ci),ci',p];
[~,p,ci] = ttest(within_VIS);
r_tab(2,:) = [mean(ci),ci',p];
[~,p,ci] = ttest(within_SOM);
r_tab(3,:) = [mean(ci),ci',p];
[~,p,ci] = ttest(between_DMN_VIS);
r_tab(4,:) = [mean(ci),ci',p];
[~,p,ci] = ttest(between_DMN_SOM);
r_tab(5,:) = [mean(ci),ci',p];
[~,p,ci] = ttest(between_VIS_SOM);
r_tab(6,:) = [mean(ci),ci',p];
[~,p,ci] = ttest(DMN_minus);
r_tab(7,:) = [mean(ci),ci',p];

r_tab = round(r_tab,2,'significant');
rownames = {'Within DMN','Within VIS','Within SOM','Between DMN-VIS',...
    'Between DMN-SOM','Between VIS-SOM','TRs: DMN-, VIS+, SOM+'};
save(fullfile(savedir,'histogramtable.mat'),'r_tab','rownames');

%%
% show each TR in DMN- on axes defined by mean activity in DMN, SOM, VIS

SOM_activity = mean(concTS(partition_data == DMN_low,netlabels_sort == SOM),2);
VIS_activity = mean(concTS(partition_data == DMN_low,netlabels_sort == VIS),2);
DMN_activity = mean(concTS(partition_data == DMN_low,netlabels_sort == DMN),2);

maxval = ceil(max([abs(SOM_activity);abs(VIS_activity);abs(DMN_activity)]));
MAE = 0.01; %alpha of points

% Fig S8e-g
f = figure;
subplot(1,3,1);
scatter(SOM_activity,VIS_activity,'MarkerEdgeAlpha',MAE);
xlim([-maxval maxval]); ylim([-maxval maxval]);
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlabel('SOM'); ylabel('VIS');
axis square

subplot(1,3,2);
scatter(DMN_activity,VIS_activity,'MarkerEdgeAlpha',MAE); 
xlim([-maxval maxval]); ylim([-maxval maxval]);
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlabel('DMN'); ylabel('VIS');
axis square

subplot(1,3,3);
scatter(DMN_activity,SOM_activity,'MarkerEdgeAlpha',MAE); 
xlim([-maxval maxval]); ylim([-maxval maxval]);
hold on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlabel('DMN'); ylabel('SOM');
axis square

f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 21 10];
f.PaperSize = [21 10];
%saveas(f,fullfile(savedir,[r_n_name,'DMN-FrameActivity.pdf']));
print(fullfile(savedir,[r_n_name,'DMN-FrameActivity.png']),'-dpng','-r300');
