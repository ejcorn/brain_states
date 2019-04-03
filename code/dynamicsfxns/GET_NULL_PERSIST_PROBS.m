function [persistenceProbability] = GET_NULL_PERSIST_PROBS(kClusterAssignments,subjInd)

% Calculate persistence probabilities for null model

kClusterAssignments = reshape(kClusterAssignments,length(kClusterAssignments),1); %convert to row vector
nobs = max(subjInd);
numClusters = length(unique(kClusterAssignments));
possible_transitions = (numClusters)*(numClusters);
persistenceProbability = zeros(nobs,numClusters);

for N = 1:nobs
    subjMask = subjInd == N;
    subjMask = kClusterAssignments(subjMask)';
    subjMask = subjMask(randperm(length(subjMask)));
    for K = 1:numClusters
        persistenceProbability(N,K) = length(strfind(subjMask,[K K])) / sum(subjMask(1:(end-1)) == K);      % divide by # of states not incl. last TR
    end
end

persistenceProbability(isnan(persistenceProbability)) = 0;	% only applies if state not present

