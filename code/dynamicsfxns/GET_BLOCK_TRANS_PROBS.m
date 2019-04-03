function [transitionProbabilityMat,transitionProbabilityVec] = GET_BLOCK_TRANS_PROBS(partition,indicator,numClusters)

% partition: Vector of integers corresponding to state time series
% indicator: binary vector indicating which states belong to which blocks
% k: number of states, i.e. 5 means assume states 1:5 exist, regardless of what exists in partition

[runs,n_runs] = GET_BLOCK_STATE_TIME_SERIES(partition,indicator);

% get number of each transition in each run of the task block
numTransitionsByBlock = cell(n_runs,1);
for j = 1:n_runs
	numTransitionsByBlock{j} = GET_NUM_TRANSITIONS(runs{j},numClusters);
end

% sum transition counts over all blocks
numTransitionsCondition = sum(cat(3,numTransitionsByBlock{:}),3);

% divide row elements by respective row sums
% such that ijth element is probability of transition from i to j given you are starting in state i and it's not the last state in the series
transitionProbabilityMat = numTransitionsCondition ./ repmat(sum(numTransitionsCondition,2),[1 numClusters]);
transitionProbabilityVec = reshape(transitionProbabilityMat',1,[]);