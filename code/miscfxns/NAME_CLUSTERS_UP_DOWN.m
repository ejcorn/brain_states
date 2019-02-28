function [clusterNames,clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(centroids)

%Provide names for clusters based on angular distance to binary Yeo
%System Vectors
%returns vector where 1 indicates a "+" state and 0 indicates a "-" state
%centroids = kClusterCentroids;

[nparc,numClusters] = size(centroids);

if nparc > 400
    load('yeo7netlabelsLaus250.mat'); network7labels = network7labels(1:nparc);
else
    load('yeo7netlabelsLaus125.mat'); network7labels = network7labels(1:nparc);
end

numNets = 7;
binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == network7labels);

YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN'};

% calculate E.D. from binary state vector to centroids

centroids_up = centroids .* (centroids > 0);
centroids_down = -1 * centroids .* (centroids < 0);     % make negative activity positive and get rid of positive activity

net7ED_Up = zeros(numClusters,numNets);
net7ED_Down = zeros(numClusters,numNets);

for K = 1:numClusters
    for B = 1:numNets
        net7ED_Up(K,B) = dot(centroids_up(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
        net7ED_Down(K,B) = dot(centroids_down(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

% get index of minimum and assign names

clusterNamesUp = cell(numClusters,1);
clusterNamesDown = cell(numClusters,1);
for K = 1:numClusters       % for up and down separately, calculate closest network
    Up_ind = find(net7ED_Up(K,:) == max(net7ED_Up(K,:)));
    Down_ind = find(net7ED_Down(K,:) == max(net7ED_Down(K,:)));    %cos(0) = 1 so need max not min (like for E.D.)
    clusterNamesUp{K} = [YeoNetNames{Up_ind},'+'];
    clusterNamesDown{K} = [YeoNetNames{Down_ind},'-'];
end

clusterNames = [clusterNamesUp,clusterNamesDown];