% script to reorder states for k = 5 so +/- states can be easily compared

a = clock; rng(a(6));
minNumClusters = 2; maxNumClusters = 11; clusterRange = minNumClusters:maxNumClusters;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(basedir,['data/TimeSeriesIndicators',name_root,'.mat']));

load(['k',num2str(numClusters),name_root,'.mat'],'clusterAssignments');

if numClusters == 5 | numClusters == 6
    clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
    if numClusters == 5
    	niceOrder = {'DMN+','DMN-','FPN+','VIS+','VIS-'};
    elseif numClusters == 6
    	niceOrder = {'DMN+','DMN-','FPN+','SOM+','VIS+','VIS-'};
    end
    shuffleIdx = zeros(numClusters,1);
    for K = 1:numClusters
    	shuffleIdx(K) = find(strcmp(clusterNames,niceOrder{K}));
    end
    
    % reorder cluster names and centroids
    clusterAssignments.(['k',num2str(numClusters)]).clusterNames = clusterNames(shuffleIdx);
    clusterAssignments.(['k',num2str(numClusters)]).clusterNamesUp = clusterAssignments.(['k',num2str(numClusters)]).clusterNamesUp(shuffleIdx);
    clusterAssignments.(['k',num2str(numClusters)]).clusterNamesDown = clusterAssignments.(['k',num2str(numClusters)]).clusterNamesDown(shuffleIdx);
    clusterAssignments.(['k',num2str(numClusters)]).bestCentroid = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid(:,shuffleIdx);
    assignments = clusterAssignments.(['k',num2str(numClusters)]).partition;
    idx = false(length(assignments),numClusters);
	
	for K = 1:numClusters
		idx(:,K) = assignments == K;
	end
	figure;subplot(1,2,1); imagesc(idx);
	for K = 1:numClusters
		assignments(idx(:,shuffleIdx(K))) = K;
	end
	idx = false(length(assignments),numClusters);
	for K = 1:numClusters
		idx(:,K) = assignments == K;
	end
	subplot(1,2,2); imagesc(idx);

	clusterAssignments.(['k',num2str(numClusters)]).partition = assignments;

	savedir = [masterdir,'/clusterAssignments'];
	cd(savedir);
	save(['k',num2str(numClusters),name_root,'.mat'],'clusterAssignments');

end


