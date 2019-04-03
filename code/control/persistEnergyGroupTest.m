addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy');
mkdir(savedir);

% empirical states
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
%% normalize states
load('data/submittedpapercentroids.mat')

StateMagnitude = sqrt(sum(kClusterCentroids.^2,1));
Xo = kClusterCentroids ./ StateMagnitude;		% normalize all states
Xf = Xo;	% start and end at same place, i.e. persistence energy

%% load structure, normalize, generate nulls
load(fullfile(savedir,['GroupRepresentativeSC_Laus',num2str(lausanneScaleBOLD),'.mat']),'A','D');

%% compute minimum control energy required to maintain each state
T = 1; % control horizon
% compute inverse of controllability gramian for each network
% normalize all networks to max eig then subtract identity to ensure max eigval = 0 so system converges on single eigenmode over time
%A = randmio_und(A,1);

%nbins = 11; nrewire = 1e5; [~,A] = fcn_preserve_degseq_lengthdist(Anorm,D,nbins,nrewire); 
%A = A + A' - eye(length(A));
WcI = GRAMIAN(A,T,true);

% compute minimum control energy required to maintain cluster centers
Epersist = MIN_CONTROL_ENERGY(A,WcI,Xo,Xf,T,true);

disp(clusterNames')
disp(Epersist)

% load persistence energy

load(fullfile(masterdir,'analyses','transitionprobabilities',['nBackCombTransitionProbabilities_k',num2str(numClusters),name_root,'.mat']),'transitionProbability')
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));	% isolate persistence probabilities from linearized transition probability matrices
restMeanPP = mean(transitionProbability(:,onDiag),1);

% correlate with persistence energy

disp('persistence probability')
corr(restMeanPP',Epersist')

load(fullfile(masterdir,'analyses','transitionprobabilities',['nBackCombDwellTime_k',num2str(numClusters),name_root,'.mat']),'DwellTimeMean')
restMeanDT = mean(DwellTimeMean,1);

disp('dwell time')
corr(restMeanDT',Epersist')