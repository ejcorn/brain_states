function [transitionProbability,transProbEndGivenStart] = GET_TRANS_PROBS(kClusterAssignments,subjInd,numClusters)

% Calculate transition probabilities from sequential cluster assignments
% Return 2D (for regression) or 3D versions of matrix (for eigen
% decompositions and plotting)

kClusterAssignments = reshape(kClusterAssignments,length(kClusterAssignments),1); %convert to row vector
nobs = max(subjInd);
if ~exist('numClusters','var')
	numClusters = length(unique(kClusterAssignments));
end
possible_transitions = (numClusters)*(numClusters);
numTransitions = zeros(nobs,numClusters,numClusters);

transProbEndGivenStart = zeros(size(numTransitions));

for N = 1:nobs
    subjMask = subjInd == N;
    subjMask = kClusterAssignments(subjMask)';
    for Kinitial = 1:numClusters
        for Kfinal = 1:numClusters
            numTransitions(N,Kinitial,Kfinal) = length(strfind(subjMask,[Kinitial Kfinal]));
        end
    end
    transProbEndGivenStart(N,:,:) = squeeze(numTransitions(N,:,:)) ./ repmat(sum(squeeze(numTransitions(N,:,:)),2),[1 numClusters]);
end

transProbEndGivenStart(isnan(transProbEndGivenStart)) = 0;

transitionProbability = reshape(permute(transProbEndGivenStart,[1 3 2]),nobs,possible_transitions);


