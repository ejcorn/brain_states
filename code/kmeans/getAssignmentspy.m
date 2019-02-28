rng('shuffle');
addpath(genpath('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode'))
addpath(genpath('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/kmeanscode'))
addpath(genpath('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/assesskmeanscode'))
addpath(genpath('/data/tesla-data/ecornblath/matlab/control_fc/NCT'));
minNumClusters = 2; maxNumClusters = 18; clusterRange = minNumClusters:maxNumClusters; numK = 1 + maxNumClusters-minNumClusters;
k_offset = minNumClusters - 1;
masterdir = ['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root];
addpath(masterdir);
load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/Demographics',name_root,'.mat']);
load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/ConcTimeSeries',name_root,'.mat'])


savedir = [masterdir,'/clusterAssignments'];
mkdir(savedir);

d = [masterdir,'/repkmeans'];
splits_by_K = zeros(nsplits,numK);
for S = 1:nsplits
  for K = 1:numK
      fname = [d,'/kmeans',num2str(S),distanceMethod,'k_',num2str(clusterRange(K)),'/kmeans',...
          num2str(S),'k_',num2str(clusterRange(K)),'rep1.mat'];
      splits_by_K(S,K) = logical(exist(fname,'file'));
  end
end

nsplits = min(sum(splits_by_K,1));

totalNumTPs = length(subjInd);


combPartitions = cell(numK,1); combSumD = cell(numK,1);
disp('start loading k means partitions');
for K = minNumClusters:maxNumClusters
    existingSplits = find(splits_by_K(:,K-k_offset));
    combPartitions{K-k_offset} = int8(zeros(totalNumTPs,nsplits));
    combSumD{K-k_offset} = zeros(1,nsplits);
        for splitInd = 1:nsplits
            S = existingSplits(splitInd);
            startInd = sum(splits_by_K(1:S,K-k_offset)) - splits_by_K(S,K-k_offset);
            load([d,'/kmeans',num2str(S),distanceMethod,'k_',num2str(K),'/kmeans',num2str(S),'k_',num2str(K),'rep1.mat']);
            combPartitions{K-k_offset}(:,splitInd) = int8(partition);
            combSumD{K-k_offset}(:,splitInd) = sum(sumd);
        end
    combPartitions{K-k_offset} = int8(multislice_pair_labeling(combPartitions{K-k_offset}));
    disp(['load partitions: K=',num2str(K)]);
end

nreps = nsplits;

cd(savedir);
save(['repkmeansPartitions',distanceMethod,name_root,'.mat'],'combPartitions','combSumD','nreps','splits_by_K','-v7.3');

concTS = dlmread(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/ConcTSCSV_',name_root,'.csv'],',');
disp('time series loaded');

disp('finding lowest MSE partition');
for numClusters = minNumClusters:maxNumClusters
    combBestClusterCentroids = zeros(nparc,numClusters);
    combBestClusterInd = find(combSumD{numClusters-k_offset} == min(combSumD{numClusters-k_offset}));
    combBestClusterInd = combBestClusterInd(1); %if multiple partitions exactly the same, just use the first one]
    
    for K = 1:numClusters
        combBestClusterCentroids(:,K) = mean(concTS(combPartitions{numClusters-k_offset}(:,combBestClusterInd) == K,:),1)';
    end
    
    clusterAssignments.(['k',num2str(numClusters)]).bestCentroid = combBestClusterCentroids;
    clusterAssignments.(['k',num2str(numClusters)]).bestClusterInd = combBestClusterInd;
    clusterAssignments.(['k',num2str(numClusters)]).partition = combPartitions{numClusters - k_offset}(:,combBestClusterInd);
    clusterAssignments.(['k',num2str(numClusters)]).sumD = combSumD{K-k_offset};
    clusterAssignments.(['k',num2str(numClusters)]).clusterNames = NAME_CLUSTERS_ANGLE(combBestClusterCentroids);
    [~,clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(combBestClusterCentroids);
    clusterAssignments.(['k',num2str(numClusters)]).clusterNamesUp = clusterNamesUp;
    clusterAssignments.(['k',num2str(numClusters)]).clusterNamesDown = clusterNamesDown;

    save(['k',num2str(numClusters),name_root,'.mat'],'clusterAssignments');
    clear clusterAssignments
    disp(['K=',num2str(K)]);
end

