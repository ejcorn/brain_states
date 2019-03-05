function [clusterNamesUp,clusterNamesDown,net7corr_Up,net7corr_Down] = NAME_CLUSTERS_CORR_UP_DOWN(centroids)

% Generate names for clusters based on correlation to binary Yeo
% This function looks at high (+) and low (-) amplitude activity separately
% and assigns a + name to the Yeo system that the regions with activity > 0 are most correlated with
% and a - name to Yeo system that regions with activity < 0 are most correlated with
% centroids: NxK matrix

[nparc,numClusters] = size(centroids);

if nparc > 400
    load('data/yeo7netlabelsLaus250.mat'); network7labels = network7labels(1:nparc);
else
    load('data/yeo7netlabelsLaus125.mat'); network7labels = network7labels(1:nparc);
end

% make a matrix where each column corresponds to a labeled Yeo system in Lausanne parcellation
% the columns are binary vectors indicated whether a region belongs to corresponding Yeo system
numNets = 7;
binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == network7labels);

YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};


% compute masks to isolate high and low amplitude regions (i.e. >0, <0)
centroids_up = centroids .* (centroids > 0);
centroids_down = centroids .* (centroids < 0); 

net7corr_Up = zeros(numClusters,numNets);
net7corr_Down = zeros(numClusters,numNets);

for K = 1:numClusters
    for B = 1:numNets
        % compute correlations with Yeo systems with low amplitude regions set to 0
        net7corr_Up(K,B) = corr(centroids_up(:,K),binaryNetVectors(:,B));
        % compute correlations with Yeo systems with high amplitude regions set to 0
        net7corr_Down(K,B) = corr(-1*centroids_down(:,K),binaryNetVectors(:,B));
    end
end




% get index of minimum and assign names

clusterNamesUp = cell(numClusters,1);
clusterNamesDown = cell(numClusters,1);
for K = 1:numClusters       % for up and down separately, calculate closest network
    Up_ind = find(net7corr_Up(K,:) == max(net7corr_Up(K,:)));
    Down_ind = find(net7corr_Down(K,:) == max(net7corr_Down(K,:)));    %cos(0) = 1 so need max not min (like for E.D.)
    clusterNamesUp{K} = [YeoNetNames{Up_ind},'+'];
    clusterNamesDown{K} = [YeoNetNames{Down_ind},'-'];
end