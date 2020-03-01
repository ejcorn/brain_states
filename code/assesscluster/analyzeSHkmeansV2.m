% makes Fig S2c-f

addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));

savedir = [masterdir,'/analyses/choosing_k/splithalves'];
mkdir(savedir);

nreps = 10;
splits_by_R = zeros(nsplits,1);
for S = 1:nsplits
  for R = 1:nreps
      fname = [savedir,'/SHkmeans',num2str(S),distanceMethod,'k_',num2str(numClusters),'/SHkmeans',...
          num2str(S),'k_',num2str(numClusters),'rep',num2str(R),'.mat'];
      splits_by_R(S,R) = logical(exist(fname,'file'));
  end
end

[splitInd,repInd] = find(splits_by_R);
totalNumTPs = length(subjInd);
nobs_train = floor(0.5*nobs);
nobs_test = nobs - nobs_train;
totalNumReps = length(splitInd);
partitions_train = zeros(nobs_train*(120 + 225),totalNumReps); % will store all training set partitions
partitions_test = zeros(nobs_test*(120 + 225),totalNumReps); % will store all training set partitions
sumD_train = zeros(totalNumReps,1);
sumD_test = zeros(totalNumReps,1);
train_idx = zeros(nobs_train,totalNumReps);
test_idx = zeros(nobs_test,totalNumReps);
disp('start loading k means partitions');
K = numClusters;
for idx = 1:totalNumReps
  S = splitInd(idx); R = repInd(idx);   % iterate through split-rep pairs to load
  load([savedir,'/SHkmeans',num2str(S),distanceMethod,'k_',num2str(K),'/SHkmeans',num2str(S),'k_',num2str(K),'rep',num2str(R),'.mat']);
  partitions_train(:,idx) = int8(partition1);
  sumD_train(idx) = sum(sumd1);
  train_idx(:,idx) = obs1;
  partitions_test(:,idx) = int8(partition2);
  sumD_test(idx) = sum(sumd2);
  test_idx(:,idx) = obs2;
end

%% centroid similarity
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));

CentroidsTestTrain = zeros(numClusters,numClusters,totalNumReps);
centroid_names_train = cell(totalNumReps,1);
centroid_names_test = cell(totalNumReps,1);
niceOrder = {'DMN+','DMN-','FPN+','VIS+','VIS-'};
for R = 1:totalNumReps
  tic
  % make train and test time series
  TS_train = concTS(ismember(subjInd,train_idx(:,R)),:);
  TS_test = concTS(ismember(subjInd,test_idx(:,R)),:);
  % calculate centroids
  centroids_train = zeros(nparc,numClusters);
  centroids_test = zeros(nparc,numClusters);
  for K = 1:numClusters
    centroids_train(:,K) = mean(TS_train(partitions_train(:,R) == K,:),1);
    centroids_test(:,K) = mean(TS_test(partitions_test(:,R) == K,:),1);
  end
  %[partitions_train(:,R),centroids_train] = REORDER_BY_NAME(partitions_train(:,R),centroids_train,niceOrder);
  %[partitions_test(:,R),centroids_test] = REORDER_BY_NAME(partitions_test(:,R),centroids_test,niceOrder);

  CentroidsTestTrain(:,:,R) = corr(centroids_train,centroids_test); %i,j is correlation of train(:,i) and train(:,j)
  disp(['Rep ',num2str(R)])
  toc
end

[maxCentroidTrainTestCorr,shuffleIdx] = max(CentroidsTestTrain,[],1);  
maxCentroidTrainTestCorr = reshape(squeeze(maxCentroidTrainTestCorr)',1,[]);
barcolors = [53 183 121; 68 1 84] / 255;

% Fig S2d
f = figure; histogram(maxCentroidTrainTestCorr,'FaceColor',barcolors(1,:),'FaceAlpha',0.5); 
ylabel('# of Centroid Pairs'); xlabel('Correlation'); xlim([-1 1]);
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
cd(savedir);
saveas(f,['SplithalvesCentroidSimilarity_k',num2str(numClusters),name_root,'.pdf'],'pdf');

%% reorder to compare trans probs
shuffleIdx = squeeze(shuffleIdx)'; % ith column of shuffleIdx indexes the train centroid most similar to the ith test centroid 

for R = 1:totalNumReps
  idx = false(length(partitions_test(:,R)),numClusters);
  % test cluster K needs to become cluster shuffleIdx(K)
  for K = 1:numClusters 
    idx(:,K) = partitions_test(:,R) == K;  % find each test cluster
  end
  for K = 1:numClusters
    partitions_test(idx(:,K),R) = shuffleIdx(R,K);  % replace with cluster # most similar to train
  end
end
partitions_test = int8(partitions_test);
partitions_train = int8(partitions_train);
save(['SplitHalvesPartitionsOrdered_k',num2str(numClusters),name_root,'.mat'],'partitions_test','partitions_train');
%% transition probability similarity, rest
numTransitions = numClusters^2;
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numTransitions); offDiag(onDiag) = [];
rest_transitions_train = zeros(totalNumReps,numTransitions);
rest_transitions_test = zeros(totalNumReps,numTransitions);

rTR = 120;
for R = 1:totalNumReps
  tmp = GET_TRANS_PROBS(partitions_train(1:(rTR*nobs_train),R),repelem(1:nobs_train,rTR),numClusters);
  rest_transitions_train(R,:) = mean(tmp,1);  
  tmp = GET_TRANS_PROBS(partitions_test(1:(rTR*nobs_test),R),repelem(1:nobs_test,rTR),numClusters);
  rest_transitions_test(R,:) = mean(tmp,1);  
  disp(['Rep ',num2str(R)])
end
test_train_TPcorr = diag(corr(rest_transitions_train(:,offDiag)',rest_transitions_test(:,offDiag)'));
test_train_PPcorr = diag(corr(rest_transitions_train(:,onDiag)',rest_transitions_test(:,onDiag)'));

% Fig S2e
f = figure; histogram(test_train_TPcorr,'FaceColor',barcolors(2,:),'FaceAlpha',0.5); 
ylabel('# of Samples'); xlabel('Correlation'); title('Transition Probability'); xlim([-1 1]);
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [4 4];
f.PaperPosition = [0 0 4 4];
cd(savedir);
saveas(f,['SplithalvesRestTPSimilarity_k',num2str(numClusters),name_root,'.pdf'],'pdf');

% Fig S2e
f = figure; histogram(test_train_PPcorr,'FaceColor',barcolors(1,:),'FaceAlpha',0.5); 
ylabel('# of Samples'); xlabel('Correlation'); title('Persistence Probability'); xlim([-1 1]);
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [4 4];
f.PaperPosition = [0 0 4 4];
cd(savedir);
saveas(f,['SplithalvesRestPPSimilarity_k',num2str(numClusters),name_root,'.pdf'],'pdf');

%% transition probability similarity, n-back

nback_transitions_train = zeros(totalNumReps,numTransitions);
nback_transitions_test = zeros(totalNumReps,numTransitions);

nTR = 225;
for R = 1:totalNumReps
  tmp = GET_TRANS_PROBS(partitions_train((1+rTR*nobs_train):end,R),repelem(1:nobs_train,nTR),numClusters);
  nback_transitions_train(R,:) = mean(tmp,1);  
  tmp = GET_TRANS_PROBS(partitions_test((1+rTR*nobs_test):end,R),repelem(1:nobs_test,nTR),numClusters);
  nback_transitions_test(R,:) = mean(tmp,1);  
  disp(['Rep ',num2str(R)])
end

test_train_TPcorr = diag(corr(nback_transitions_train(:,offDiag)',nback_transitions_test(:,offDiag)'));
test_train_PPcorr = diag(corr(nback_transitions_train(:,onDiag)',nback_transitions_test(:,onDiag)'));

% Fig S2f
f = figure; histogram(test_train_TPcorr,'FaceColor',barcolors(2,:),'FaceAlpha',0.5); 
ylabel('# of Samples'); xlabel('Correlation'); title('Transition Probability'); xlim([-1 1]);
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [4 4];
f.PaperPosition = [0 0 4 4];
cd(savedir);
saveas(f,['SplithalvesnbackTPSimilarity_k',num2str(numClusters),name_root,'.pdf'],'pdf');

% Fig S2f
f = figure; histogram(test_train_PPcorr,'FaceColor',barcolors(1,:),'FaceAlpha',0.5); 
ylabel('# of Samples'); xlabel('Correlation'); title('Persistence Probability'); xlim([-1 1]);
prettifyEJC; 
f.PaperUnits = 'centimeters';
f.PaperSize = [4 4];
f.PaperPosition = [0 0 4 4];
cd(savedir);
saveas(f,['SplithalvesnbackPPSimilarity_k',num2str(numClusters),name_root,'.pdf'],'pdf');