function [numTransMat,numTransVec] = GET_NUM_TRANSITIONS_NO_PERSIST(partition,numClusters)

% partition: Vector of integers corresponding to state time series
% numClusters: number of states, i.e. 5 means assume states 1:5 exist, regardless of what exists in partition
% calculate transition probabilities while excluding state persistence, i.e. independent of autocorrelation

if ~exist('numClusters','var')
	numClusters = length(unique(kClusterAssignments));
end

partition = reshape(partition,1,[]);		% make partition row vector for strfind to work
partition = [partition(find(diff(partition) ~= 0)),partition(end)]; % exclude temporal persistence to allow transitions to be at any time lag

numTransMat = zeros(numClusters);

for Kinitial = 1:numClusters
    for Kfinal = 1:numClusters
    	% calculate number of times one cluster precedes another
        numTransMat(Kinitial,Kfinal) = length(strfind(partition,[Kinitial Kfinal]));
    end
end

% flatten by row
transitionProbability = reshape(numTransMat',1,[]);