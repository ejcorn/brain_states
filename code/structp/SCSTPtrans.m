% IMPORTANT: this version 'double counts' connections between regions that belong to both states. i.e it counts them once as coming from node i, and once
% as coming from node j. Were there asymmetry, this would be the correct way to perform this analysis.

% There is also a second block of code that allows you to entirely exclude connections between regions belonging to both states

addpath(genpath('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode'));
addpath(genpath('/data/tesla-data/ecornblath/matlab/BCT/'));

masterdir = ['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root];
load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/Demographics',name_root,'.mat']);

load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat']);
kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;

load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/VolNormSC',num2str(lausanneScaleBOLD),'.mat']);
savedir = [masterdir,'/analyses/structrans'];
mkdir(savedir);

possible_transitions = numClusters^2;
corrlog = logical(triu(true(nparc,nparc)) .* ~eye(nparc));

centroidKeyNodes = zscore(kClusterCentroids) > thrsh;

interStateSC = zeros(nobs,numClusters,numClusters);
interStateSTP = zeros(nobs,numClusters,numClusters);
interStateComm = zeros(nobs,numClusters,numClusters);tic
overlap = zeros(numClusters);
for N = 1:nobs
    % put nan's on diagonal so you don't include it in averaging
    % this will only matter if there is overlap in active nodes between states
    A = SCvolnorm{N};  
    STP = distance_bin(A>0) + diag(nan(nparc,1));
    D = diag(strengths_und(A))^-0.5;
    Comm = expm(D*A*D) + diag(nan(nparc,1));
    A = A + diag(nan(nparc,1));
    for K1 = 1:numClusters
        T1 = find(centroidKeyNodes(:,K1));
        for K2 = 1:numClusters
            T2 = find(centroidKeyNodes(:,K2));             
            %{
            % 'double count' connections between regions in both states
            interStateSC(N,K1,K2) = nanmean(nanmean(A(T1,T2)));
            interStateSTP(N,K1,K2) = nanmean(nanmean(STP(T1,T2)));
            interStateComm(N,K1,K2) = nanmean(nanmean(Comm(T1,T2)));        
            %}
            % exclude all connections between regions in both states
            
            A_s = A; % make new matrices with NaNs on all connections between regions in both states
            Comm_s = Comm;
            STP_s = STP;
            overlap(K1,K2) = sum(ismember(T1,T2));
            A_s(T1(ismember(T1,T2)),T2(ismember(T2,T1))) = NaN;
            A_s(T2(ismember(T2,T1)),T1(ismember(T1,T2))) = NaN;
            STP_s(T1(ismember(T1,T2)),T2(ismember(T2,T1))) = NaN;
            Comm_s(T1(ismember(T1,T2)),T2(ismember(T2,T1))) = NaN;

            interStateSC(N,K1,K2) = nanmean(reshape(A_s(T1,T2),1,[]));
            interStateSTP(N,K1,K2) = nanmean(reshape(STP_s(T1,T2),1,[]));
            interStateComm(N,K1,K2) = nanmean(reshape(Comm_s(T1,T2),1,[]));
            
            %{
            % only look at connections between regions in both states
            interStateSC(N,K1,K2) = nanmean(nanmean(A(T1(ismember(T1,T2)),T2(ismember(T2,T1)))));
            interStateSTP(N,K1,K2) = nanmean(nanmean(STP(T1(ismember(T1,T2)),T2(ismember(T2,T1)))));
            interStateComm(N,K1,K2) = nanmean(nanmean(Comm(T1(ismember(T1,T2)),T2(ismember(T2,T1)))));       
            %}

        end
    end
    disp(['Subject = ',num2str(N)]);
end
toc
interStateSC = reshape(permute(interStateSC,[1 3 2]),nobs,possible_transitions);
interStateSTP = reshape(permute(interStateSTP,[1 3 2]),nobs,possible_transitions);
interStateComm = reshape(permute(interStateComm,[1 3 2]),nobs,possible_transitions);


cd(savedir);
save(['structrans_k',num2str(numClusters),'thresh',num2str(thrsh),name_root,'.mat'],'interStateSC','interStateSTP','interStateComm','overlap');