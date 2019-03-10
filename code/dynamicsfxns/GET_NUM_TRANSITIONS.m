function [numTransMat,numTransVec] = GET_NUM_TRANSITIONS(partition,numClusters)

% partition: Vector of integers corresponding to state time series
% numClusters: number of states, i.e. 5 means assume states 1:5 exist, regardless of what exists in partition

if ~exist('numClusters','var')
	numClusters = length(unique(kClusterAssignments));
end

partition = reshape(partition,1,[]);		% make partition row vector for strfind to work
numTransMat = zeros(numClusters);

for Kinitial = 1:numClusters
    for Kfinal = 1:numClusters
    	% calculate number of times one cluster precedes another
        numTransMat(Kinitial,Kfinal) = length(strfind(partition,[Kinitial Kfinal]));
    end
end

% flatten by row
transitionProbability = reshape(numTransMat',1,[]);