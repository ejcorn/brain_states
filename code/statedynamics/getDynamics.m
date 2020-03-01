addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(basedir,['data/TimeSeriesIndicators',name_root,'.mat']));

if scan == 'R'
    scanlab = {'Rest'}; numTRs = 120;
elseif scan == 'N'
    scanlab = {'nBack'}; numTRs = 225;
    scanInd = scanInd * 0;  %because scanInd is all 1's 
elseif scan == 'C'
    scanlab = {'RestComb','nBackComb'}; numTRs = [120 225]; % index rest num TRs (120) and nback to loop through
end

savedir = [masterdir,'/analyses/transitionprobabilities'];
mkdir(savedir);
load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

for i = 1:numel(scanlab)
    
    %% calculate transition prob

    [transitionProbability,transitionProbabilityMats,numTransitions] = GET_TRANS_PROBS(kClusterAssignments(scanInd == (i-1)),subjInd(scanInd == (i-1)));    
    save(fullfile(savedir,[scanlab{i},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']),'transitionProbabilityMats')
    save(fullfile(savedir,[scanlab{i},'TransitionProbabilities_k',num2str(numClusters),name_root,'.mat']),'transitionProbability')
    save(fullfile(savedir,[scanlab{i},'NumTransitions_k',num2str(numClusters),name_root,'.mat']),'numTransitions')

    [transitionProbability,transitionProbabilityMats] = GET_TRANS_PROBS_NO_PERSIST(kClusterAssignments(scanInd == (i-1)),subjInd(scanInd == (i-1)));    
    save(fullfile(savedir,[scanlab{i},'TransitionProbabilitiesNoPersist_k',num2str(numClusters),name_root,'.mat']),'transitionProbability','transitionProbabilityMats');

    %% calculate fractional occupancy, i.e. percentage of time spent in each state, regardless of order
    
    FractionalOccupancy = zeros(nobs,numClusters);

    for N = 1:nobs
        for K = 1:numClusters
            FractionalOccupancy(N,K) = sum(and(and(kClusterAssignments == K,subjInd' == N),scanInd == (i-1))) / numTRs(i);
        end
    end

    save(fullfile(savedir,[scanlab{i},'FractionalOccupancy_k',num2str(numClusters),name_root,'.mat']),'FractionalOccupancy')

    %% calculate dwell time, i.e. average number of subsequent TRs that each state lasts for
    % store both mean and median, because dwell time may not be normally distributed
    TR = 3;     % PNC TR length
    DwellTimeMean = zeros(nobs,numClusters);
    DwellTimeMedian = zeros(nobs,numClusters);
    RunRate = zeros(nobs,numClusters);
    for N = 1:nobs
        [dt_mean,dt_median,~,~,n_runs] = CALC_DWELL_TIME(kClusterAssignments(subjInd' == N & scanInd == (i-1)),numClusters);
        DwellTimeMean(N,:) = dt_mean*TR;        % store dwell time in seconds
        DwellTimeMedian(N,:) = dt_median*TR;        % store dwell time in seconds
        % store rate of appearance of runs, i.e. how many DMN runs appear in 1 minute
        % first calculate runs/TR,then runs/sec, then runs/min
        RunRate(N,:) = 60*(1/TR)*(n_runs/numTRs(i));   
    end

    save(fullfile(savedir,[scanlab{i},'DwellTime_k',num2str(numClusters),name_root,'.mat']),'DwellTimeMean','DwellTimeMedian','RunRate')

end

