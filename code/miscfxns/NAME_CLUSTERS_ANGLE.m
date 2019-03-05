function [clusterNames,reorderClusters,clusterNamesSort,net7angle] = NAME_CLUSTERS_ANGLE(centroids)

%Provide names for clusters based on angular distance to binary Yeo
%System Vectors
%returns vector where 1 indicates a "+" state and 0 indicates a "-" state
%centroids = kClusterCentroids;

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

% calculate cosine of angle between binary state vector and centroids

net7angle = zeros(numClusters,numNets*2);

for K = 1:numClusters
    for B = 1:(numNets*2)
        net7angle(K,B) = dot(centroids(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

% get index of minimum and assign names

clusterNamesInit = cell(numClusters,1);
plusminus = true(numClusters,1);
for K = 1:numClusters
    ind = find(net7angle(K,:) == max(net7angle(K,:)));    %cos(0) = 1 so need max not min (like for E.D.)
    clusterNamesInit{K} = YeoNetNames{ind};
    if ind > numNets
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
