rng(split);
nreps = 10;
addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
savedir = fullfile(masterdir,'analyses','choosing_k','splithalves');
mkdir(savedir); cd(savedir);

concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));
if zdim > 0
    concTS = zscore(concTS,[],zdim);
end

disp('data loaded');

disp(['K = ',num2str(numClusters),'Split = ',num2str(split)]);

savedir = ['SHkmeans',num2str(split),distanceMethod,'k_',num2str(numClusters)];
mkdir(savedir);
%%
cd(savedir); 

disp('start k-means');
for R = 1:nreps
	obs1 = randperm(nobs,floor(0.5*nobs));
    [partition1,~,sumd1] = kmeans(concTS(ismember(subjInd,obs1),:),numClusters,'Distance',distanceMethod);
    partition1 = int8(partition1);

    obs2 = find(~ismember(1:nobs,obs1));	% cluster the other half
    [partition2,~,sumd2] = kmeans(concTS(ismember(subjInd,obs2),:),numClusters,'Distance',distanceMethod);
    partition2 = int8(partition2); 
    save(['SHkmeans',num2str(split),'k_',num2str(numClusters),'rep',num2str(R),'.mat'],'partition1','sumd1','obs1','partition2','sumd2','obs2');
    clear partition1; clear partition2;
    clear sumd1; clear sumd2;
    disp(['K-means ',num2str(R)]);
end

disp('complete');
