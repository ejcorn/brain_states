function [pvals_twotail] = PERM_TEST(A,B,nperms)

    % test null hypothesis that group average transition matrix elements in A do not differ from group average of B
    % by swapping 50% of the membership between A and B and recomputing difference between A and B

    % INPUT:
    % A and B: NxKxK stack of KxK transition matrices for N subjects and two conditions, A and B
    % nperms: number of permutations, pretty fast so 10000 is good to get a precise p-value
        
    % OUTPUT:
    % pvals_twotail: equal tail bootstrap p-value for A ~= B, 

    % NOTES:
    % using nanmean in the event that a transition is not observed in a subject
    % we handle this by ignoring the transition rather than setting it to 0
    % this does not occur at all if every subject has every state in the state time series
    % used to calculate transition probs, which is the case for rest and n-back scans at k=5.
    % however, within n-back blocks, there are some subjects that don't have every state
    
    disp('start permutation testing')
    [nobs, numClusters, numClusters] = size(A);
    permutationIndex = rand(nobs,nperms) > 0.5; % preallocate uniform random swap indices
    expectedA_gt_B = zeros(nperms,numClusters,numClusters); % difference under null hypothesis 
    for P = 1:nperms
        shuffleIndex = find(permutationIndex(:,P));
        BTmp = B; ATmp = A;
        
        % Randomly switch B and n-back entire matrices within each subject
        BTmp(shuffleIndex,:,:) = A(shuffleIndex,:,:);
        ATmp(shuffleIndex,:,:) = B(shuffleIndex,:,:);

        % Compute amount by which mean B exceeds mean nback 
        expectedA_gt_B(P,:,:) = squeeze(nanmean(ATmp,1) - nanmean(BTmp,1));
        disp(['Perm ',num2str(P)])
    end

    %% compare null distribution to observed 
    disp('compare null to observed')

    % Compute actual difference between mean B and mean A
    observedA_gt_B = repmat(nanmean(A,1) - nanmean(B,1),[nperms 1 1]);

    % Compute equal tail p-value for (A - B) > null
    % this should be plotted on the actual A-B
    % 
    A_gt_null = squeeze(sum(observedA_gt_B > expectedA_gt_B,1)) / nperms;  % number of times observed difference greater than null difference
    A_lt_null = squeeze(sum(observedA_gt_B < expectedA_gt_B,1)) / nperms;  % number of times observed difference less than null difference
    pvals_twotail = 2*min(A_gt_null,A_lt_null);
    