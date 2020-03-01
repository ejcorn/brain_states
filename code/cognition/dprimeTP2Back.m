addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','cognition','transprobs');
mkdir(savedir);

%%
load(['data/Demographics',name_root,'.mat']);
load(fullfile(savedir,['TwoBackTPNoPersistDprime_k',num2str(numClusters),'.mat'])); % load regression output
load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

f=figure;

subplot(1,3,1);
imagesc(bmat); h =colorbar;
ylabel(h,'\beta_{TP}');
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
[y,x] = find(pmat < 0.05); % already bonferroni corrected in R
text(x-.12,y+.12,'*','Color','w');
ylabel('Current State'); xlabel('Next State');
title('2-back Performance');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
save(fullfile(savedir,'Fig4d__TPvs2backDPrime.mat'),'bmat','pmat');
saveas(f,fullfile(savedir,['TwoBackTPNoPersistDprimeMatlab_k',num2str(numClusters),'.pdf']));
