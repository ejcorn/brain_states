addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));

scanlab = {'RestComb','nBackComb'}; 

savedir = [masterdir,'/analyses/transitionprobabilities'];

rng('shuffle');

%% load transition probabilities
disp('loading transition probabilities')

load([masterdir,'/analyses/transitionprobabilities/',scanlab{1},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
restTransitionProbabilityMats = transitionProbabilityMats; clear transitionProbabilityMats
load([masterdir,'/analyses/transitionprobabilities/',scanlab{2},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
nBackTransitionProbabilityMats = transitionProbabilityMats; clear transitionProbabilityMats

load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

%% permute rest and n-back within subjects
disp('start permutation testing')
nperms = 10000;
permutationIndex = rand(nobs,nperms) > 0.5;
expectedRest_gt_nBack = zeros(nperms,numClusters,numClusters);
expectednBack_gt_Rest = zeros(nperms,numClusters,numClusters);
for P = 1:nperms
    shuffleIndex = find(permutationIndex(:,P));
    restTmp = restTransitionProbabilityMats; nBackTmp = nBackTransitionProbabilityMats;
    
    % Randomly switch rest and n-back entire matrices within each subject
    restTmp(shuffleIndex,:,:) = nBackTransitionProbabilityMats(shuffleIndex,:,:);
    nBackTmp(shuffleIndex,:,:) = restTransitionProbabilityMats(shuffleIndex,:,:);
    
    % Count # of times rest exceeds nback in each permutation
    %expectedRest_gt_nBack(P,:,:) = squeeze(sum(restTmp > nBackTmp,1));
    %expectednBack_gt_Rest(P,:,:) = squeeze(sum(restTmp < nBackTmp,1));    

    % Compute amount by which mean rest exceeds mean nback 
    expectednBack_gt_Rest(P,:,:) = squeeze(mean(nBackTmp,1) - mean(restTmp,1));

end

%% compare null distribution to observed 
disp('compare null to observed')
% Count # of times rest exceeds nback in actual data
%observedRest_gt_nBack = repmat(sum(restTransitionProbabilityMats > nBackTransitionProbabilityMats,1),[nperms 1 1]);
%observednBack_gt_Rest = repmat(sum(restTransitionProbabilityMats < nBackTransitionProbabilityMats,1),[nperms 1 1]);

% Compute actual difference between mean rest and mean nback
observednBack_gt_Rest = repmat(mean(nBackTransitionProbabilityMats,1) - mean(restTransitionProbabilityMats,1),[nperms 1 1]);

% Count % of times observed exceeds each null value
%Rest_gt_nBack_Pval = squeeze(sum(observedRest_gt_nBack > expectedRest_gt_nBack,1)) / nperms;
%nBack_gt_Rest_Pval = squeeze(sum(observednBack_gt_Rest > expectednBack_gt_Rest,1)) / nperms;

% Compute percentage of times observed difference exceeded null
nBack_gt_Rest_Pval = squeeze(sum(observednBack_gt_Rest > expectednBack_gt_Rest,1)) / nperms;


%% Plot 
disp('plot')

cd(savedir);

grpAvgRest = squeeze(mean(restTransitionProbabilityMats,1));
grpAvgRest = grpAvgRest .* ~eye(numClusters);
grpAvgnBack = squeeze(mean(nBackTransitionProbabilityMats,1)); %restGrpAvg(logical(eye(numClusters))) = 0;
grpAvgnBack = grpAvgnBack .* ~eye(numClusters);
maxval = max(max([grpAvgRest,grpAvgnBack]))

% add asterisks for significance relative to null
lt = 0.05/(numClusters^2 - numClusters); ut = 1-lt;   % Bonferroni corrected thresholds -- transitions only not persistence
load(['randnull/',scanlab{1},'Null_ShuffledStatesD_k',num2str(numClusters),'.mat'])
[ygt,xgt] = find((probExceedNull > ut).*~eye(numClusters));
[ylt,xlt] = find((probExceedNull < lt).*~eye(numClusters));

f=figure;
subplot(1,3,1);
imagesc(grpAvgRest); 
text(xgt,ygt,'+','Color','w');
text(xlt,ylt,'-','Color','w');
caxis([0 maxval]); colorbar
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title('Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

clear probExceedNull
load(['randnull/',scanlab{2},'Null_ShuffledStatesD_k',num2str(numClusters),'.mat'])
[ygt,xgt] = find((probExceedNull > ut).*~eye(numClusters));
[ylt,xlt] = find((probExceedNull < lt).*~eye(numClusters));
cd(savedir);
subplot(1,3,2); 
imagesc(grpAvgnBack); 
text(xgt,ygt,'+','Color','w');
text(xlt,ylt,'-','Color','w');
caxis([0 maxval]); colorbar
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title('n-back');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

subplot(1,3,3);
imagesc((nBack_gt_Rest_Pval.*~eye(numClusters))); colormap('plasma');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('State at t'); xlabel('State at t + 1');
lt = 0.05/(numClusters^2); ut = 1-lt;       % bonferroni correction
[y,x] = find(nBack_gt_Rest_Pval.*~eye(numClusters) > ut | nBack_gt_Rest_Pval.*~eye(numClusters) < lt);
text(x-.12,y+.12,'*','Color','w');
%h=colorbar; ylabel(h,'-log_{10}(p)'); caxis([0 log10(nperms)]); h.Ticks = [0 ut lt log10(nperms)]; h.TickLabels = [0 round(ut,2,'significant') round(lt,2,'significant') log10(nperms)];
h = colorbar; ylabel(h,'p(n > R)'); caxis([lt ut]); h.Ticks = [lt 0.5 ut]; h.TickLabels = [round(lt,2,'significant') 0.5 round(ut,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('n-back > Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

saveas(f,['RestvsnBack_nonpar_k',num2str(numClusters),'.pdf']);
[y,x] = find(diag(nBack_gt_Rest_Pval)' > ut | diag(nBack_gt_Rest_Pval)' < lt);
f=figure;
imagesc(diag(nBack_gt_Rest_Pval)');
xticks([]); yticks([]);
text(x-.12,y+.12,'*','Color','w');
caxis([lt ut]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];
saveas(f,['PersistRestvsnBack_nonpar_k',num2str(numClusters),'.pdf']);