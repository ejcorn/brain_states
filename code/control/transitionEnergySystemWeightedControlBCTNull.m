addpaths;
%% load matrix, define real initial and final conditions
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
savedir = fullfile(basedir,'results/',...
    name_root,'analyses','control_energy','model_task_inputs',['BCTNull_c',num2str(c),'T',num2str(T)]);

mkdir(savedir);

% load yeo labels
load(fullfile(basedir,'data',['yeo7netlabelsLaus',num2str(lausanneScaleBOLD),'.mat']));
SystemLabels = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN','SUB'};
InputSystems = find(ismember(SystemLabels,SystemLabel)); % get indices of systems whose control inputs should be weighted
InputVector = ismember(network7labels(1:nparc),InputSystems);
disp(['Control inputs weighted 2:1 towards ',num2str(sum(InputVector)),' nodes in the ', SystemLabel,' system.'])
B = repmat(InputVector,[1 nparc]) .*eye(nparc); % construct input matrix allowing input only into selected regions in InputVector

% option to bias input towards input vector but allow input into all regions
if strcmp('Weighted',ControlLabel) 
   B = B + eye(nparc); 
elseif strcmp('Weighted10fold',ControlLabel) % now test differences in arbitrary scaling between input and background controllers
   B = B + 0.1*eye(nparc); 
elseif strcmp('Weighted1_10',ControlLabel)
   B = 9*B + eye(nparc); 
elseif strcmp('Weighted1_5',ControlLabel)
   B = 4*B + eye(nparc); 
elseif strcmp('Weighted2_4',ControlLabel) % keep fold difference in input weight constant
   B = 2*(B + eye(nparc));
elseif strcmp('Weighted0.5_1.5',ControlLabel) % keep absolute difference in input weight constant
   B = (B + eye(nparc)) - 0.5;
end

load(['results/',name_root,'/analyses/centroids/OverallClusterCentroids_k',num2str(numClusters),name_root,'.mat']);
Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = REPELEM_SUBSTITUTE(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix

x0 = kClusterCentroids(:,Xo_ind);
xf = kClusterCentroids(:,Xf_ind);

num_transitions = numClusters^2;

% load group SC
load(fullfile(basedir,['results/',name_root,'/analyses/control_energy/GroupRepresentativeSC_FA_Laus',num2str(lausanneScaleBOLD),'.mat']),'A');

%% compute control energy for state transitions while controlling from yeo systems

n_per_split = ceil(nperms/nsplits); % calculate subjects per job based on # of jobs (nsplits) and # of subjects (nobs)
perm_range = (1+n_per_split*(split-1)):(n_per_split*split);

E_full = zeros(num_transitions,1);
DistanceToTarget = zeros(num_transitions,1);
NumericalError = zeros(num_transitions,1);
for P = perm_range
  disp(['Perm ',num2str(P)])
  A_BCTnull = randmio_und(A,1);
  PreservedDegreeSequence = isequal(degrees_und(A),degrees_und(A_BCTnull)); % check if model preserved degree sequence
  A_BCTnullnorm = NORMALIZE(A_BCTnull,c); % stabilize null matrix  
  for transition = 1:num_transitions    
       [x, u, NumericalError(transition)] = MIN_ENG_CONT(A_BCTnullnorm, T, B, x0(:,transition), xf(:,transition), 0);       
       DistanceToTarget(transition) = sum((xf(:,transition)-x(end,:)').^2);
       E_full(transition) = sum(sum(u.^2))*T/1000;
  end
  save(fullfile(savedir,['BCTNullPerm',num2str(P),'_FA_TransitionEnergies',SystemLabel,ControlLabel,'Control_c',...
  num2str(c),'T',num2str(T),'_k',num2str(numClusters),'.mat']),'E_full','DistanceToTarget','NumericalError','PreservedDegreeSequence');
end