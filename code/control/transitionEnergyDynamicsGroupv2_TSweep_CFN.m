addpaths;
%% load matrix, define real initial and final conditions
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
% load group representative SC
load(fullfile(basedir,['results/',name_root,'/analyses/control_energy/GroupRepresentativeSC_FA_Laus',num2str(lausanneScaleBOLD),'.mat']));
savedir = fullfile(masterdir,'analyses',...
 'control_energy','Tsweep',['c',num2str(c),'_k',num2str(numClusters)]); % for cluster
mkdir(savedir);

Anorm = NORMALIZE(A,c);

load(['results/',name_root,'/analyses/centroids/OverallClusterCentroids_k',num2str(numClusters),name_root,'.mat']);
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix

x0 = kClusterCentroids(:,Xo_ind);
xf = kClusterCentroids(:,Xf_ind);
%% load null states
% make null states
nreps = 20; nsplits = 50;
nperms = nreps*nsplits;
num_transitions = numClusters^2;

%% compute control energy for state transitions, real and null
T_rng = [0.001:0.5:10]; nT = length(T_rng);
E_full = zeros(num_transitions,nT);
E_null_full = zeros(num_transitions,nperms,nT);
xf_to_natural_endpoint = zeros(num_transitions,nT); % distance of natural trajectory to endpoint

for T_i = 1:nT
    T = T_rng(T_i);
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    E_full(:,T_i) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false);    
    xf_to_natural_endpoint(:,T_i) = sum((xf - expm(Anorm*T)*x0).^2,1); % distance to natural trajectory of system
end

%% make a new A with the same eigenvectors as structure but skewed eigenvalues:
gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2))
skewedgaussian = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x)

%% compute control energy for state transitions in null networks

nperms = 500;
E_BCTnull = zeros(num_transitions,nperms,nT);
%E_URnull = zeros(num_transitions,nperms,nT); % uniform random networks
%E_LSnull = zeros(num_transitions,nperms,nT); % networks with same eigenvectors as brain but very different eigenvalues
%E_VVnull = zeros(num_transitions,nperms,nT); % networks with random eigenvectors and different distribution of eigenvalues

%E_BCTnull_nullstates = zeros(numClusters^2,nperms);
PerservedDegreeSequence = zeros(nperms,1);
for P = 1:nperms
    disp(['BCT Perm ',num2str(P)]);
    Arand = randmio_und(A,1);
    PerservedDegreeSequence(P) = isequal(degrees_und(Arand),degrees_und(A)); % double check that degree sequence was preserved
    Arandnorm = NORMALIZE(Arand,c); % stabilize null matrices
    %{
    % other null models:
    ArandUR = rand(nparc); 
    %ArandURnorm = NORMALIZE(ArandUR + ArandUR',c);
    % make new matrices with same eigenvectors as A but new eigenvalues
    [U,D,V] = eig(A);
    D(~~eye(size(A,1))) = skewedgaussian(linspace(0,0.5+(randn(1)/10),nparc),3); % generate matrices with variably skewed eigenvalues
    A_lambdaskew = NORMALIZE(U*D*inv(U),c);
    
    % make new matrices with random eigenvectors AND a different eigenvalue distribution    
    [U,D,V] = eig(ArandUR + ArandUR');
    D(~~eye(size(A,1))) = skewedgaussian(linspace(0,0.5,nparc),3); % generate matrices with variably skewed eigenvalues
    A_valvecskew = NORMALIZE(U*D*inv(U),c);
    %}
    for T_i = 1:nT
        T = T_rng(T_i);
        WcIrand = GRAMIAN_FAST(Arandnorm, T); % compute gramian inverse for control horizon T
        E_BCTnull(:,P,T_i) = MIN_CONTROL_ENERGY(Arandnorm, WcIrand, x0, xf, T,false);    
        %WcIrand = GRAMIAN_FAST(ArandURnorm, T); % compute gramian inverse for control horizon T
        %E_URnull(:,P,T_i) = MIN_CONTROL_ENERGY(ArandURnorm, WcIrand, x0, xf, T,false);    
        %WcIrand = GRAMIAN_FAST(A_lambdaskew, T);
        %E_LSnull(:,P,T_i) = MIN_CONTROL_ENERGY(A_lambdaskew, WcIrand, x0, xf, T,false);    
        %WcIrand = GRAMIAN_FAST(A_valvecskew, T);
        %E_VVnull(:,P,T_i) = MIN_CONTROL_ENERGY(A_valvecskew, WcIrand, x0, xf, T,false);
    end
end

%% save energies

save(fullfile(savedir,['GroupAverageTransitionEnergies_k',num2str(numClusters),'.mat']),...
    'E_full','E_BCTnull','xf_to_natural_endpoint','T_rng','PerservedDegreeSequence'); %,'E_URnull','E_LSnull','E_VVnull');
