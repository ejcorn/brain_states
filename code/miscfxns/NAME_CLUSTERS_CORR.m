function [clusterNames,reorderClusters,clusterNamesSort,net7corr,sysCorrs] = NAME_CLUSTERS_CORR(centroids)

% Generate names for clusters based on the largest absolute correlation
% between the centroid and either positive or negative binary indicators
% of Yeo systems.
% centroids: NxK matrix
% net7corr: KxY matrix containing correlation values between centroids and each Yeo system, positive or negative
% sysCorrs: Kx1 vector containing correlation values with closeest system for each centroid.

[nparc,numClusters] = size(centroids);

if nparc > 400
    load('data/yeo7netlabelsLaus250.mat'); network7labels = network7labels(1:nparc);
else
    load('data/yeo7netlabelsLaus125.mat'); network7labels = network7labels(1:nparc);
end

numNets = 7;
% make a matrix where each column corresponds to a labeled Yeo system in Lausanne parcellation
% the columns are binary vectors indicated whether a region belongs to corresponding Yeo system

binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == network7labels);
% then duplicate this matrix, multiply by -1 and horizontally concatenate to
% provide separate names for when systems are low amplitude

binaryNetVectors = [binaryNetVectors, -1*binaryNetVectors];

YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-'};

% calculate correlation between each binary state vector and each centroid

net7corr = zeros(numClusters,numNets*2);

for K = 1:numClusters
    for B = 1:(numNets*2)
        net7corr(K,B) = corr(centroids(:,K),binaryNetVectors(:,B));
    end
end

% get index of minimum and assign names

clusterNamesInit = cell(numClusters,1);
plusminus = true(numClusters,1);
sysCorrs = zeros(numClusters,1);                 
for K = 1:numClusters
    ind = find(net7corr(K,:) == max(net7corr(K,:)));    % find maximum correlation for each cluster
    clusterNamesInit{K} = YeoNetNames{ind};             % give that cluster the corresponding name
    sysCorrs(K) = net7corr(K,ind);      % save correlation value with most similar binary vector
    if ind > 7
        plusminus(K) = false;
    end
end

clusterNames = cellstr(clusterNamesInit);

%sort by name then plus-minus
%
[clusterNamesSort,I] = sort(clusterNamesInit);
[~,I2] = sort(plusminus(I));
clusterNamesSort = clusterNamesSort(I2);
reorderClusters = I(I2);
%}
