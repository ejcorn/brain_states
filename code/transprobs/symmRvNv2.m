addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));

if scan ~= 'C'
	return
end
scanlab = {'RestComb','nBackComb'}; numTRs = [120 225]; % index rest num TRs (120) and nback to loop through

savedir = [masterdir,'/analyses/transitionprobabilities/symmetry'];
mkdir(savedir);
cd(savedir);

a = clock;
rng(a(6));

%% load transition probabilities
disp('loading transition probabilities')

load([masterdir,'/analyses/transitionprobabilities/',scanlab{1},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
transitionProbabilityMatsRest= transitionProbabilityMats; clear transitionProbabilityMats
load([masterdir,'/analyses/transitionprobabilities/',scanlab{2},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
transitionProbabilityMatsnBack = transitionProbabilityMats; clear transitionProbabilityMats


RestNull_UniformGivenSubjectDiag = zeros(nobs,numClusters,numClusters);
nBackNull_UniformGivenSubjectDiag = zeros(nobs,numClusters,numClusters);
for N = 1:nobs

    PersistenceProbs = squeeze(transitionProbabilityMatsRest(N,:,:)) .* eye(numClusters);
    % make matrix with 1 - P(persist in k) / numClusters as off diagonals in row k
    RestNull_UniformGivenSubjectDiag(N,:,:) = (ones(numClusters) * (eye(numClusters)-PersistenceProbs)/(numClusters-1))' .* ~eye(numClusters) + PersistenceProbs;
    
    PersistenceProbs = squeeze(transitionProbabilityMatsnBack(N,:,:)) .* eye(numClusters);
    % make matrix with 1 - P(persist in k) / numClusters as off diagonals in row k
    nBackNull_UniformGivenSubjectDiag(N,:,:) = (ones(numClusters) * (eye(numClusters)-PersistenceProbs)/(numClusters-1))' .* ~eye(numClusters) + PersistenceProbs;

end


%[~,transitionProbabilityMatsNull] = GET_TRANS_PROBS(tmpAssignments,subjInd(scanInd == 0));

restSymmetryScore = zeros(nobs,1);
nBackSymmetryScore = zeros(nobs,1);
restNullSymmetryScore = zeros(nobs,1);
nBackNullSymmetryScore = zeros(nobs,1);

for N = 1:nobs
	restSymmetryScore(N) = SYMMETRY_SCORE(transitionProbabilityMatsRest(N,:,:));
	nBackSymmetryScore(N) = SYMMETRY_SCORE(transitionProbabilityMatsnBack(N,:,:));
	restNullSymmetryScore(N) = SYMMETRY_SCORE(RestNull_UniformGivenSubjectDiag(N,:,:));
	nBackNullSymmetryScore(N) = SYMMETRY_SCORE(nBackNull_UniformGivenSubjectDiag(N,:,:));
end

save(['RvNSymmetryScorev2_k',num2str(numClusters),'.mat'],'restSymmetryScore','nBackSymmetryScore','restNullSymmetryScore','nBackNullSymmetryScore');



