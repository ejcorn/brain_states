addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
savedir = [masterdir,'/analyses/transitionprobabilities'];
scanlab = {'RestComb','nBackComb'};

load([masterdir,'/analyses/transitionprobabilities/',scanlab{1},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
restTransitionProbabilityMats = transitionProbabilityMats; clear transitionProbabilityMats
load([masterdir,'/analyses/transitionprobabilities/',scanlab{2},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
nBackTransitionProbabilityMats = transitionProbabilityMats; clear transitionProbabilityMats
load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;

rPP = diag(squeeze(mean(restTransitionProbabilityMats,1)))';
nPP = diag(squeeze(mean(nBackTransitionProbabilityMats,1)))';
maxval = max([rPP,nPP])
nperms = 2500; 

%% rest
disp('starting to calculate shuffled persistence probabilities Rest');
tic
tmpAssignments = kClusterAssignments(scanInd == 0);	% rest
subjIndtmp = subjInd(scanInd == 0);
probExceedNull = zeros(nperms,numClusters);
parfor P = 1:nperms
    tic
    nullPersistenceProbs = GET_NULL_PERSIST_PROBS(tmpAssignments,subjIndtmp);
    toc
    tmp = rPP > mean(nullPersistenceProbs,1);        
    probExceedNull(P,:) = tmp;
    disp(['Perm ',num2str(P)])
end
disp('done shuffling transition probabilities');
toc
probExceedNull = mean(probExceedNull,1);	% p-value

lt = 0.025/(2*numClusters^2); ut = 1-lt;   % two-tailed Bonferroni corrected thresholds over whole matrix
[ygt,xgt] = find((probExceedNull > ut));
[ylt,xlt] = find((probExceedNull < lt));

f=figure;
imagesc(rPP);
xticks([]); yticks([]);
text(xgt,ygt,'+','Color','w');
text(xlt,ylt,'-','Color','w');
caxis([0 maxval]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];
cd(savedir);
saveas(f,['RestPersistenceProbs_k',num2str(numClusters),'.pdf']);

% for source data file
if strcmp(name_root,'ScanCLaus250Z0final') && numClusters == 5
    save(['Fig4a__RestPersistProbs_k',num2str(numClusters),'.mat'],'rPP','probExceedNull');
elseif strcmp(name_root,'ScanCLaus250Z0final') && numClusters == 4
    save(['FigS13c__RestPersistProbs_k',num2str(numClusters),'.mat'],'rPP','probExceedNull');
elseif strcmp(name_root,'ScanCLaus250Z0final') && numClusters == 6
    save(['FigS14c__RestPersistProbs_k',num2str(numClusters),'.mat'],'rPP','probExceedNull');
elseif strcmp(name_root,'ScanCLaus250Z0cosinefinal') && numClusters == 5
    save(['FigS12c__RestPersistProbs_k',num2str(numClusters),'.mat'],'rPP','probExceedNull');
elseif strcmp(name_root,'ScanCLaus125Z0final') && numClusters == 5
    save(['FigS15c__RestPersistProbs_k',num2str(numClusters),'.mat'],'rPP','probExceedNull');
end

%% n-back
disp('starting to calculate shuffled persistence probabilities n-back');
tic

tmpAssignments = kClusterAssignments(scanInd == 1);	% rest
subjIndtmp = subjInd(scanInd == 1);
probExceedNull = zeros(nperms,numClusters);
parfor P = 1:nperms
    tic
    nullPersistenceProbs = GET_NULL_PERSIST_PROBS(tmpAssignments,subjIndtmp);
    toc
    tmp = nPP > mean(nullPersistenceProbs,1);        
    probExceedNull(P,:) = tmp;
    disp(['Perm ',num2str(P)])
end
disp('done shuffling transition probabilities');
toc
probExceedNull = mean(probExceedNull,1);	% p-value

lt = 0.025/(2*numClusters^2); ut = 1-lt;   % two-tailed Bonferroni corrected thresholds over whole matrix
[ygt,xgt] = find((probExceedNull > ut));
[ylt,xlt] = find((probExceedNull < lt));

f=figure;
imagesc(nPP);
xticks([]); yticks([]);
text(xgt,ygt,'+','Color','w');
text(xlt,ylt,'-','Color','w');
caxis([0 maxval]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];
cd(savedir);
saveas(f,['nBackPersistenceProbs_k',num2str(numClusters),'.pdf']);

% for source data file

if strcmp(name_root,'ScanCLaus250Z0final') && numClusters == 5
    save(['Fig4b__nBackPersistProbs_k',num2str(numClusters),'.mat'],'nPP','probExceedNull');
elseif strcmp(name_root,'ScanCLaus250Z0final') && numClusters == 4
    save(['FigS13d__nBackPersistProbs_k',num2str(numClusters),'.mat'],'nPP','probExceedNull');
elseif strcmp(name_root,'ScanCLaus250Z0final') && numClusters == 6
    save(['FigS14d__nBackPersistProbs_k',num2str(numClusters),'.mat'],'nPP','probExceedNull');
elseif strcmp(name_root,'ScanCLaus250Z0cosinefinal') && numClusters == 5
    save(['FigS12d__nBackPersistProbs_k',num2str(numClusters),'.mat'],'nPP','probExceedNull');
elseif strcmp(name_root,'ScanCLaus125Z0final') && numClusters == 5
    save(['FigS15d__nBackPersistProbs_k',num2str(numClusters),'.mat'],'nPP','probExceedNull');
end
