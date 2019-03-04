addpaths;
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','hcpLR');
nparc = 462;
rTR = round(405); nTR = 405;

[hcprTP,HCPresttransprobs] = GET_TRANS_PROBS(partitionPNCorder(scanInd == 0),HCPsubjInd(scanInd == 0));
[hcpnTP,HCPnBacktransprobs] = GET_TRANS_PROBS(partitionPNCorder(scanInd == 1),HCPsubjInd(scanInd == 1));
hcprTP = mean(hcprTP,1);
hcpnTP = mean(hcpnTP,1);

% compare HCP and PNC 
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = [];

load([masterdir,'/analyses/transitionprobabilities/RestCombTransitionProbabilities_k',num2str(numClusters),name_root,'.mat'])
pncrTP = mean(transitionProbability,1);
load([masterdir,'/analyses/transitionprobabilities/nBackCombTransitionProbabilities_k',num2str(numClusters),name_root,'.mat'])
pncnTP = mean(transitionProbability,1);

% correlate TPs and PPs between HCP and PNC
restTPsim = round(corr(pncrTP(offDiag)',hcprTP(offDiag)'),2,'significant');
restPPsim = round(corr(pncrTP(onDiag)',hcprTP(onDiag)'),2,'significant');
nbackTPsim = round(corr(pncnTP(offDiag)',hcpnTP(offDiag)'),2,'significant');
nbackPPsim = round(corr(pncnTP(onDiag)',hcpnTP(onDiag)'),2,'significant');

% fold difference in TP and PP
mean(pncrTP(offDiag) ./ hcprTP(offDiag))
mean(hcprTP(onDiag) ./ pncrTP(onDiag))
mean(pncnTP(offDiag) ./ hcpnTP(offDiag))
mean(hcpnTP(onDiag) ./ pncnTP(onDiag))

%% test rest vs. n-back transitions

restTransitionProbabilityMats = HCPresttransprobs;
nBackTransitionProbabilityMats = HCPnBacktransprobs;

%% permute rest and n-back within subjects
disp('start permutation testing')
nperms = 10000;
permutationIndex = rand(nsubjs,nperms) > 0.5;
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
observednBack_gt_Rest = repmat(mean(nBackTransitionProbabilityMats,1) - mean(restTransitionProbabilityMats,1),[nperms 1 1]);
nBack_gt_Rest_Pval = squeeze(sum(observednBack_gt_Rest > expectednBack_gt_Rest,1)) / nperms;

% plot
HCPgrpAvgRest = squeeze(mean(HCPresttransprobs,1)) .* ~eye(numClusters);
HCPgrpAvgnBack = squeeze(mean(HCPnBacktransprobs,1)) .* ~eye(numClusters);
maxVal = max(max([HCPgrpAvgRest,HCPgrpAvgnBack]));

f = figure;

subplot(1,3,1);
imagesc(HCPgrpAvgRest);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title({'Rest',['r = ',num2str(restTPsim),' with PNC']});
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,3,2);
imagesc(HCPgrpAvgnBack);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);

ylabel('State at t'); xlabel('State at t + 1');
title({'n-back',['r = ',num2str(nbackTPsim),' with PNC']});
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

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

saveas(f,fullfile(savedir,['HCP_XH_RestnBackTransProbs_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.pdf']),'pdf');

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
saveas(f,fullfile(savedir,['HCP_XH_RestvsnBackPersistProbs_k',num2str(numClusters),'.pdf']));

%% Plot persistence probabilities
rPP = diag(squeeze(mean(restTransitionProbabilityMats,1)))';
nPP = diag(squeeze(mean(nBackTransitionProbabilityMats,1)))';
maxval = max([rPP,nPP])

f=figure;
imagesc(rPP);
xticks([]); yticks([]);
caxis([0 maxval]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];

saveas(f,fullfile(savedir,['HCP_XH_RestPersistenceProbs_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),'.pdf']),'pdf');

f=figure;
imagesc(nPP);
xticks([]); yticks([]);
caxis([0 maxval]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];

saveas(f,fullfile(savedir,['HCP_XH_nBackPersistenceProbs_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),'.pdf']),'pdf');

