function [dwellTimesByCluster] = GET_BLOCK_DWELL_TIME(partition,indicator,numClusters)

% partition: Vector of integers corresponding to state time series
% indicator: binary vector indicating which states belong to which blocks
% k: number of states, i.e. 5 means assume states 1:5 exist, regardless of what exists in partition

[runs,n_runs] = GET_BLOCK_STATE_TIME_SERIES(partition,indicator);

% get dwell times in each run of the task block
dwellTimesByBlock = cell(n_runs,1);

for j = 1:n_runs
	[~,~,~,dwellTimesByBlock{j}] = CALC_DWELL_TIME(runs{j},numClusters);
end

% get dwell times by cluster for each block
for k = 1:numClusters
	dwellTimesByCluster{k} = cell(n_runs,1);
	for j = 1:n_runs
		dwellTimesByCluster{k}{j} = dwellTimesByBlock{j}{k};
	end
	dwellTimesByCluster{k} = cat(2,dwellTimesByCluster{k}{:});
end

% get mean dwell time across all runs
dwellTimesByCluster = cellfun(@mean,dwellTimesByCluster);