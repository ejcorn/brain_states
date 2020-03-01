function [transitionProbabilityMat,transitionProbabilityVec,numTransitionsMat,numTransitionsVec] = GET_BETWEEN_BLOCK_TRANS_PROBS(partition,indicator,numClusters)

% partition: Vector of integers corresponding to state time series
% indicator: vector with unique integers indicating which TRs belong to which blocks
% k: number of states, i.e. 5 means assume states 1:5 exist, regardless of what exists in partition

blockSwitches = find(diff(indicator) ~= 0); % indexes the last TR of each block, so that index and that index + 1 highlight TRs on either side of transition boundaries
numTransitionsMat = zeros(numClusters);
for BS = blockSwitches'
	state_i = partition(BS); % state at end of a block
	state_j = partition(BS+1); % state at beginning of next block
	numTransitionsMat(state_i,state_j) = numTransitionsMat(state_i,state_j) + 1; % add to count of transitions
end

numTransitionsVec = reshape(numTransitionsMat',1,[]);

% divide row elements by respective row sums
% such that ijth element is probability of transition from i to j given you are starting in state i and it's not the last state in the series
transitionProbabilityMat = numTransitionsMat ./ repmat(sum(numTransitionsMat,2),[1 numClusters]);
transitionProbabilityVec = reshape(transitionProbabilityMat',1,[]);