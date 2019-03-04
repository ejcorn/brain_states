a=clock;
rng(a(6));

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
savedir = fullfile(masterdir,'repkmeans');
mkdir(savedir); cd(savedir);
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));


if zdim > 0
	concTS = zscore(concTS,[],zdim);
end

disp('data loaded');

disp(['K = ',num2str(numClusters),'Split = ',num2str(split)]);

savedir = ['kmeans',num2str(split),distanceMethod,'k_',num2str(numClusters)];
mkdir(savedir);
%%
cd(savedir);

disp('start k-means');
for R = 1:nreps
    [partition,~,sumd] = kmeans(concTS,numClusters,'Distance',distanceMethod);
    save(['kmeans',num2str(split),'k_',num2str(numClusters),'rep',num2str(R),'.mat'],'partition','sumd');
    clear partition
    clear sumd
    disp(['K-means ',num2str(R)]);
end

disp('complete');
