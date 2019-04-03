addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
if scan == 'R'
    scanlab = {'Rest'}; numTRs = 120;
    scanttl = {'Rest'};
elseif scan == 'N'
    scanlab = {'nBack'}; numTRs = 225;
    scanInd = scanInd * 0;  %because scanInd is all 1's 
    scanttl = {'n-back'};
elseif scan == 'C'
    scanlab = {'RestComb','nBackComb'}; numTRs = [120 225]; % index rest num TRs (120) and nback to loop through
    scanttl = {'Rest','n-back'};
end
savedir = fullfile(masterdir,'/analyses/transitionprobabilities/randnull');
mkdir(savedir);

a = clock; rng(a(6));

for i = 1:numel(scanlab)

    load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
    kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;
    kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
    clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

    %load([masterdir,'/analyses/transitionprobabilities/',scanlab{i},'TransitionProbabilities_k',num2str(numClusters),name_root,'.mat']);

    tmpAssignments = kClusterAssignments(scanInd == (i-1))';    %analyze rest and n-back one at a time
    possible_transitions = numClusters^2;

    subjIndtmp = subjInd(scanInd == (i-1));
    transitionProbability = GET_TRANS_PROBS_NO_PERSIST(tmpAssignments,subjIndtmp);
    probExceedNull = zeros(nperms,possible_transitions);
    grpAvg = mean(transitionProbability,1);
    disp(['starting to calculate shuffled transition probabilities ',scanlab{i}]);
    tic
    parfor P = 1:nperms
        tic
        nullTransitionProbability = GET_NULL_TRANS_PROBS(tmpAssignments,subjIndtmp);
        toc
        tmp = grpAvg > mean(nullTransitionProbability,1);        
        probExceedNull(P,:) = tmp;
        disp(['Perm ',num2str(P)])
    end
    disp('done shuffling transition probabilities');
    toc

    probExceedNull = reshape(mean(probExceedNull,1),numClusters,numClusters)';
    
    f=figure;
    imagesc((probExceedNull)); colormap('plasma');
    xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
    yticks(1:numClusters); yticklabels(clusterNames); axis square
    ylabel('State at t'); xlabel('State at t + 1');
    lt = 0.05/(numClusters^2 - numClusters); ut = 1-lt;   % Bonferroni corrected thresholds -- transitions only not persistence
    [ygt,xgt] = find((probExceedNull > ut).*~eye(numClusters));
    [ylt,xlt] = find((probExceedNull < lt).*~eye(numClusters));
    text(xgt,ygt,'+','Color','w');
    text(xlt,ylt,'-','Color','w');
    %lt = -log10(lt); ut = -log10(ut);
    %h=colorbar; ylabel(h,'-log_{10}(p)'); caxis([0 log10(nperms)]); h.Ticks = [0 ut lt log10(nperms)]; h.TickLabels = [0 round(ut,2,'significant') round(lt,2,'significant') log10(nperms)];
    h = colorbar; ylabel(h,'p > null'); caxis([lt ut]); h.Ticks = [lt 0.5 ut]; h.TickLabels = [round(lt,2,'significant') 0.5 round(ut,2,'significant')];
    COLOR_TICK_LABELS(true,true,numClusters);
    title([scanttl{i},' vs. Null']);
    set(gca,'FontSize',8);
    f.PaperUnits = 'inches';
    f.PaperSize = [8/3 2];
    f.PaperPosition = [0 0 8/3 2];
    saveas(f,fullfile(savedir,[scanlab{i},'TPvsNull_k',num2str(numClusters),'.pdf']));
    save(fullfile(savedir,[scanlab{i},'Null_ShuffledStatesD_k',num2str(numClusters),'.mat']),'probExceedNull');

end
