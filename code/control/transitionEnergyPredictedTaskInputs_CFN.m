%%
% this code runs on CfN
% test whether using predicted task inputs explains n-back state
% transitions better than equal input into all regions
addpath(genpath(fullfile(basedir,'code')));

%% load matrix, define real initial and final conditions
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
% load group representative SC
load(fullfile(basedir,['results/',name_root,'/analyses/control_energy/GroupRepresentativeSC_FA_Laus250.mat']));
addpath(genpath('~/Dropbox/Cornblath_Bassett_Projects/code/BCT/'));
savedir = fullfile(basedir,'results/',...
    name_root,'analyses','control_energy','model_task_inputs');
%savedir = fullfile(masterdir,'analyses',...
% 'control_energy','Tsweep',['c',num2str(c),'_k',num2str(numClusters)]); % for cluster
mkdir(savedir);

% load yeo labels
load(fullfile(basedir,'data',['yeo7netlabelsLaus',num2str(lausanneScaleBOLD),'.mat']));

InputVector = ismember(network7labels(1:nparc),InputSystems); % only give input into these regions
B = InputVector .*eye(nparc); % construct input matrix allowing input only into selected regions in InputVector

% option to bias input towards input vector but allow input into all regions
ControlLabel = 'Weighted'; % setting this as fixed here but could theoretically pass in to try a range
if strcmp('Weighted',ControlLabel) 
   B = B + eye(nparc); 
elseif strcmp('Weighted10fold',ControlLabel) % now test differences in arbitrary scaling between input and background controllers
   B = B + 0.1*eye(nparc); 
elseif strcmp('Weighted1_10',ControlLabel)
   B = 9*B + eye(nparc); 
elseif strcmp('Weighted1_5',ControlLabel)
   B = 4*B + eye(nparc); 
elseif strcmp('Weighted4_8',ControlLabel)
   B = 4*(B + eye(nparc)); 
elseif strcmp('Weighted2_4',ControlLabel) % keep fold difference in input weight constant
   B = 2*(B + eye(nparc));
elseif strcmp('Weighted0.5_1.5',ControlLabel) % keep absolute difference in input weight constant
   B = (B + eye(nparc)) - 0.5;
elseif strcmp('Betas',ControlLabel) % use 2b-0b GLM betas
    X = load(fullfile(basedir,'results/',...
    name_root,'analyses','control_energy',['LausanneBetas',num2str(lausanneScaleBOLD),'.mat']))
    B = diag(X.lausannebetas(1:nparc));
end

%A = rand(nparc); A = A + A';
Anorm = NORMALIZE(A,c);

load(['results/',name_root,'/analyses/centroids/OverallClusterCentroids_k',num2str(numClusters),name_root,'.mat']);
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix

x0 = kClusterCentroids(:,Xo_ind);
xf = kClusterCentroids(:,Xf_ind);

num_transitions = numClusters^2;
%% compute control energy for state transitions while controlling from yeo systems

T_rng = [0.001:0.5:10]; nT = length(T_rng);
E_full = zeros(num_transitions,nT);
E_identity = zeros(num_transitions,nT); % control from all nodes
DistanceToTarget = zeros(num_transitions,nT);
NumericalError = zeros(num_transitions,nT);
%Inputs = zeros(1001,nparc,num_transitions,nT);
for T_i = 1:nT
    T = T_rng(T_i); disp(['T #',num2str(T_i)]);
    for transition = 1:num_transitions    
         [x, u, NumericalError(transition,T_i)] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
         disp(['numerical error = ',num2str(NumericalError(transition,T_i))]);
         DistanceToTarget(transition,T_i) = sum((xf(:,transition)-x(end,:)').^2);
         disp(['dist to target =',num2str(DistanceToTarget(transition,T_i))]);
         E_full(transition,T_i) = sum(sum(u.^2))*T/1001;
         Inputs(:,:,transition,T_i) = u;
    end
    WcI = GRAMIAN_FAST(Anorm,T);
    E_identity(:,T_i) = MIN_CONTROL_ENERGY(Anorm,WcI,x0,xf,T,false);
end

%% compute control energy for state transitions in null networks
%{
nperms = 100;
E_BCTnull = zeros(num_transitions,nperms,nT);
%E_BCTnull_nullstates = zeros(numClusters^2,nperms);
PerservedDegreeSequence = zeros(nperms,1);
for P = 1:nperms
    disp(['BCT Perm ',num2str(P)]);
    Arand = randmio_und(A,100);
    PerservedDegreeSequence(P) = isequal(degrees_und(Arand),degrees_und(A));
    Arandnorm = NORMALIZE(Arand,c);
    for T_i = 1:nT
        T = T_rng(T_i);
        WcIrand = GRAMIAN_FAST(Arandnorm, T); % compute gramian inverse for control horizon T
        E_BCTnull(:,P,T_i) = MIN_CONTROL_ENERGY(Arandnorm, WcIrand, x0, xf, T,false);    
    end
    %E_BCTnull_nullstates(:,P) = MIN_CONTROL_ENERGY(Arandnorm, WcIrand, x0_rand, xf_rand, T,false);
end
%}
%% save energies

save(fullfile(savedir,['TransitionEnergies_c',num2str(c),SystemLabel,ControlLabel,'Control_k',num2str(numClusters),'.mat']),...
    'E_full','E_identity','T_rng');%,'DistanceToTarget','NumericalError','Inputs');
