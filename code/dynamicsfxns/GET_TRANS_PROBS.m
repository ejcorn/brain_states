function [transitionProbability,transitionProbabilityMatrices,numTransitions] = GET_TRANS_PROBS(partition,subjInd,numClusters)

% Calculate transition probabilities in a sequence of states
% partition: integer vector, sequential cluster assignments 
% subjInd: integer vector, subject index for partition
% numClusters: number of states (indexed 1:numClusters)
% Return 2D (transitionProbability) or 3D (transitionProbabilityMatrices)

partition = reshape(partition,length(partition),1); %convert to row vector
nobs = max(subjInd);
if ~exist('numClusters','var')
	numClusters = length(unique(partition));
end
possible_transitions = (numClusters)*(numClusters);
numTransitions = zeros(nobs,numClusters,numClusters);

transitionProbabilityMatrices = zeros(size(numTransitions));

for N = 1:nobs
    subjMask = subjInd == N;
    subjMask = partition(subjMask)';
    for Kinitial = 1:numClusters
        for Kfinal = 1:numClusters
            numTransitions(N,Kinitial,Kfinal) = length(strfind(subjMask,[Kinitial Kfinal]));
        end
    end
    transitionProbabilityMatrices(N,:,:) = squeeze(numTransitions(N,:,:)) ./ repmat(sum(squeeze(numTransitions(N,:,:)),2),[1 numClusters]);
end

transitionProbabilityMatrices(isnan(transitionProbabilityMatrices)) = 0;

transitionProbability = reshape(permute(transitionProbabilityMatrices,[1 3 2]),nobs,possible_transitions);

numTransitions = reshape(permute(numTransitions,[1 3 2]),nobs,possible_transitions);