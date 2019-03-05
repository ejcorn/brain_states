addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));

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

    [transitionProbability,transitionProbabilityMats] = GET_TRANS_PROBS(kClusterAssignments(scanInd == (i-1)),subjInd(scanInd == (i-1)));
    possible_transitions = numClusters^2;    
    save(fullfile(savedir,[scanlab{i},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']),'transitionProbabilityMats')
    save(fullfile(savedir,[scanlab{i},'TransitionProbabilities_k',num2str(numClusters),name_root,'.mat']),'transitionProbability')

    %% calculate fractional occupancy, i.e. percentage of time spent in each state, regardless of order
    
    FractionalOccupancy = zeros(nobs,numClusters);

    for N = 1:nobs
        for K = 1:numClusters
            FractionalOccupancy(N,K) = sum(and(and(kClusterAssignments == K,subjInd' == N),scanInd == (i-1))) / numTRs(i);
        end
    end

    save(fullfile(savedir,[scanlab{i},'FractionalOccupancy_k',num2str(numClusters),name_root,'.mat']),'FractionalOccupancy')

    %$ calculate dwell time, i.e. average number of subsequent TRs that each state lasts for
    % store both mean and median, because dwell time may not be normally distributed
    TR = 3;     % PNC TR length
    DwellTimeMean = zeros(nobs,numClusters);
    DwellTimeMedian = zeros(nobs,numClusters);
    for N = 1:nobs
        [dt_mean,dt_median] = CALC_DWELL_TIME(kClusterAssignments(subjInd' == N & scanInd == (i-1)),numClusters);
        DwellTimeMean(N,:) = dt*TR;        % store dwell time in seconds
        DwellTimeMedian(N,:) = dt*TR;        % store dwell time in seconds
    end

    save(fullfile(savedir,[scanlab{i},'DwellTime_k',num2str(numClusters),name_root,'.mat']),'DwellTimeMean','DwellTimeMedian')

end

