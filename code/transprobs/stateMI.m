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

savedir = fullfile(masterdir,'/analyses/transitionprobabilities');
mkdir(savedir);

load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

for i = 1:numel(scanlab)   

    NGAP = 1; maxGap = 1;
    % calculate MI for each subject between states at time t and time t + k for k = 0:maxGap
    subjectMI = zeros(nobs,1);
    subjectMITransition = zeros(nobs,1);
    for N = 1:nobs    	
        mask = subjInd' == N & scanInd == (i-1);
        a = kClusterAssignments(mask);
        
        t_plus_k = [zeros(NGAP,1);a];
        % truncate state series by maxGap so all MI calcs are done on same length vectors
        t = a((maxGap):end); t_plus_k = t_plus_k((maxGap):(end-NGAP));
        subjectMI(N) = MUTUAL_INFO(t,t_plus_k);
        % shuffle t, reconstruct t_plus_k
        %

        a = [a(find(diff(a) ~= 0));a(end)]; % only assess transitions
        t_plus_k = [zeros(NGAP,1);a];
        % truncate state series by maxGap so all MI calcs are done on same length vectors
        t = a((maxGap):end); t_plus_k = t_plus_k((maxGap):(end-NGAP));
        subjectMITransition(N) = MUTUAL_INFO(t,t_plus_k);

        disp(['Subject ',num2str(N)])         
    end
    save(fullfile(savedir,[scanlab{i},'MutualInfo_k',num2str(numClusters),name_root,'.mat']),'subjectMI','subjectMITransition');
end