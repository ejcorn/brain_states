addpath(genpath('/data/tesla-data/ecornblath/matlab/control_fc/pipeline/analysiscode'));
addpath(genpath('/data/tesla-data/ecornblath/matlab/BCT/'));

masterdir = ['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root];
load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/Demographics',name_root,'.mat']);

load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/clusterTransitions_',name_root,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat']);
kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;

load(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/VolNormSCRickNull',num2str(lausanneScaleBOLD),'.mat']);
savedir = [masterdir,'/analyses/structrans'];
mkdir(savedir);

possible_transitions = numClusters^2;
corrlog = triu(true(nparc,nparc)); corrlog(logical(eye(nparc))) = 0; corrlog = logical(corrlog);

centroidKeyNodes = zscore(kClusterCentroids) > thrsh;

interStateSC = zeros(nobs,numClusters,numClusters);
interStateSTP = zeros(nobs,numClusters,numClusters);
interStateComm = zeros(nobs,numClusters,numClusters);

for N = 1:nobs
    tic
    A = SCvolnormNULL{N};    
    A1 = A(corrlog);
    A2 = distance_bin(A>0); A2 = A2(corrlog);
    Di = diag(strengths_und(A))^-0.5;
    A3 = expm(Di*A*Di); A3 = A3(corrlog);
    for K1 = 1:numClusters
        T1 = find(centroidKeyNodes(:,K1));
        
        for K2 = 1:numClusters
            T2 = find(centroidKeyNodes(:,K2)); 
            B = zeros(nparc); B(T1,T2) = 1; B(T2,T1) = 1; B = logical(B(corrlog));
            SCtmp = A1(B); %SCtmp = SCtmp(SCtmp ~= 0);
            interStateSC(N,K1,K2) = mean(SCtmp);
            STPtmp = A2(B); %STPtmp = STPtmp(STPtmp ~= 0);
            interStateSTP(N,K1,K2) = mean(STPtmp);
            Commtmp = A3(B); %Commtmp = Commtmp(Commtmp ~= 0);
            interStateComm(N,K1,K2) = mean(Commtmp);     
        end
    end
    disp(['Subject = ',num2str(N)]);
    toc
end

interStateSC = reshape(permute(interStateSC,[1 3 2]),nobs,possible_transitions);
interStateSTP = reshape(permute(interStateSTP,[1 3 2]),nobs,possible_transitions);
interStateComm = reshape(permute(interStateComm,[1 3 2]),nobs,possible_transitions);


cd(savedir);
save(['NULLstructrans_k',num2str(numClusters),'thresh',num2str(thrsh),name_root,'.mat'],'interStateSC','interStateSTP','interStateComm');