%{
% use this start to run locally
basedir=pwd; name_root = 'ScanCLaus250Z0final'; numClusters = 5; nsplits = 25; nperms = 1000;
addpath(genpath(pwd)); 
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
basedir=pwd;
masterdir = fullfile(basedir,'results',name_root);
%savedir = fullfile(masterdir,'analyses','control_energy');
savedir = fullfile(masterdir,'control_energy','sweep');
mkdir(savedir);
%}
%%
%
addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy','sweep');
mkdir(savedir);
%}
n_per_split = ceil(nperms/nsplits);	% calculate subjects per job based on # of jobs (nsplits) and # of subjects (nobs)

%% load persistence energies for control horizon and normalization sweep

load(fullfile(savedir,['SweepHorizonNormalizationValueRanges.mat']),'T_rng','c_rng');
n_horizons = length(T_rng);
n_normalizations = length(c_rng);

sweepPersistenceEnergy = zeros(n_horizons,n_normalizations,numClusters,nperms);
for split = 1:nsplits	% only store arrays at size of number of perms you process at once, then later reload in order
	perm_range = (1+n_per_split*(split-1)):(n_per_split*split);
	splitSweep = load(fullfile(savedir,['SweepHorizonNormalizationPersistenceEnergy_k',num2str(numClusters),'Split',num2str(split),'.mat']),'sweepPersistenceEnergy');
	sweepPersistenceEnergy(:,:,:,perm_range) = splitSweep.sweepPersistenceEnergy;
end

%% plot mean energy and variance

load(fullfile(masterdir,'analyses','transitionprobabilities',...
    ['OverallClusterCentroids_k',num2str(numClusters),name_root,'.mat']),'clusterNames');
% mean energy
f = figure; 
avgEnergy = mean(sweepPersistenceEnergy,4);
for K = 1:numClusters
    subplot(1,numClusters,K);    
    imagesc(avgEnergy(:,:,K)); axis square; colorbar;
    yticks(1:4:n_horizons); yticklabels(T_rng(1:4:n_horizons)); 
    xticks(1:4:n_normalizations); xticklabels(c_rng(1:4:n_normalizations));
    xlabel('c'); ylabel('T'); title(clusterNames{K});
    set(gca,'FontSize',8);
end
f.PaperUnits = 'centimeters';
f.PaperSize = [21 4];
f.PaperPosition = [0 0 21 4];
saveas(f,fullfile(savedir,['MeanPersistEnergySweep_k',num2str(numClusters),'.pdf']));

% variance of energy
f = figure; 
avgEnergy = std(sweepPersistenceEnergy,[],4);
for K = 1:numClusters
    subplot(1,numClusters,K);    
    imagesc(avgEnergy(:,:,K)); axis square; colorbar;
    yticks(1:4:n_horizons); yticklabels(T_rng(1:4:n_horizons)); 
    xticks(1:4:n_normalizations); xticklabels(c_rng(1:4:n_normalizations));
    xlabel('c'); ylabel('T'); title(clusterNames{K});
    set(gca,'FontSize',8);
end
f.PaperUnits = 'centimeters';
f.PaperSize = [21 4];
f.PaperPosition = [0 0 21 4];
saveas(f,fullfile(savedir,['VariancePersistEnergySweep_k',num2str(numClusters),'.pdf']));

% coefficient of variation of energy
f = figure; 
avgEnergy = std(sweepPersistenceEnergy,[],4) ./ mean(sweepPersistenceEnergy,4);
for K = 1:numClusters
    subplot(1,numClusters,K);    
    imagesc(avgEnergy(:,:,K)); axis square; colorbar;
    yticks(1:4:n_horizons); yticklabels(T_rng(1:4:n_horizons)); 
    xticks(1:4:n_normalizations); xticklabels(c_rng(1:4:n_normalizations));
    xlabel('c'); ylabel('T'); title(clusterNames{K});
    set(gca,'FontSize',8);
end
f.PaperUnits = 'centimeters';
f.PaperSize = [21 4];
f.PaperPosition = [0 0 21 4];
saveas(f,fullfile(savedir,['CVPersistEnergySweep_k',num2str(numClusters),'.pdf']));

%% correlate with dwell time

% load dwell time
rDT = load(fullfile(masterdir,'analyses','transitionprobabilities',...
    ['RestCombDwellTime_k',num2str(numClusters),name_root,'.mat']));
nDT = load(fullfile(masterdir,'analyses','transitionprobabilities',...
    ['nBackCombDwellTime_k',num2str(numClusters),name_root,'.mat']));

% load subject indices for each bootstrapped sample to bootstrap dwell time
bootsamps = zeros(nperms,nobs);
for P = 1:nperms
    disp(['Perm ',num2str(P)])
    load(fullfile(datadir,'BootstrapGroupSC',['BootstrapGroupRepresentativeSC',num2str(lausanneScaleBOLD),'Perm',num2str(P),'.mat']),'bootsamp');
    bootsamps(P,:) = bootsamp;  % save this to bootstrap persistence probabilities/dwell times as well
end

sweepCorrRest = zeros(n_horizons,n_normalizations,nperms);
sweepCorrnBack = zeros(n_horizons,n_normalizations,nperms);

for P = 1:nperms
    % get bootstrapped dwell time for each sample
    rDT_i = mean(rDT.DwellTimeMean(bootsamps(P,:),:),1)';
    nDT_i = mean(nDT.DwellTimeMean(bootsamps(P,:),:),1)';
    for T_i = 1:n_horizons		
		for c_i = 1:n_normalizations
            Pe_i = squeeze(sweepPersistenceEnergy(T_i,c_i,:,P));
            sweepCorrRest(T_i,c_i,P) = corr(Pe_i,rDT_i);
            sweepCorrnBack(T_i,c_i,P) = corr(Pe_i,nDT_i);
        end
    end            
end

%% plot correlations
% one-tailed tests to see if rest correlation is negative
% and n-back correlation is positive

% correct over entire sweep
bonf_thresh = 0.05/(n_horizons*n_normalizations);
pvalRest = sum(sweepCorrRest > 0,3) < bonf_thresh;
pvalnBack = sum(sweepCorrnBack < 0,3) < bonf_thresh;

f = figure; 
subplot(1,2,1);
avgRestCorr = mean(sweepCorrRest,3); % average over bootstraps
imagesc(avgRestCorr); axis square; colorbar; 
lim = avgRestCorr(find(abs(avgRestCorr) == max(max(abs(avgRestCorr)))));
colormap('plasma'); caxis(sort([0 lim]));
[yp,xp] = find(pvalRest);
text(xp,yp,'*','Color','w','FontSize',6);
yticks(1:4:n_horizons); yticklabels(T_rng(1:4:n_horizons)); 
xticks(1:4:n_normalizations); xticklabels(c_rng(1:4:n_normalizations));
xlabel('c'); ylabel('T');
title('Rest');
set(gca,'FontSize',8);

subplot(1,2,2);
avgnBackCorr = mean(sweepCorrnBack,3); % average over bootstraps
imagesc(avgnBackCorr); axis square; colorbar;
lim = avgnBackCorr(find(abs(avgnBackCorr) == max(max(abs(avgnBackCorr)))));
colormap('plasma'); caxis(sort([0 lim]));
[yp,xp] = find(pvalnBack);
text(xp,yp,'*','Color','w','FontSize',6);
yticks(1:4:n_horizons); yticklabels(T_rng(1:4:n_horizons)); 
xticks(1:4:n_normalizations); xticklabels(c_rng(1:4:n_normalizations));
xlabel('c'); ylabel('T');
title('n-back');
set(gca,'FontSize',8);

f.PaperUnits = 'centimeters';
f.PaperSize = [9 6];
f.PaperPosition = [0 0 9 6];
saveas(f,fullfile(savedir,['DT_PESweepCorrMean_k',num2str(numClusters),'.pdf']));

%% save correlations for plotting

sweepCorrRest = reshape(sweepCorrRest,[],1);
sweepCorrnBack = reshape(sweepCorrnBack,[],1);
save(fullfile(savedir,['SweepDTPECorrs_k',num2str(numClusters),'.mat']),'sweepCorrRest','sweepCorrnBack');
