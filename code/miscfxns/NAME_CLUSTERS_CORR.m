function [clusterNames,reorderClusters,clusterNamesSort] = NAME_CLUSTERS_CORR(centroids,sclfactor)

%Provide names for clusters based on angular distance to binary Yeo
%System Vectors
%returns vector where 1 indicates a "+" state and 0 indicates a "-" state
%centroids = kClusterCentroids;

if ~exist('sclfactor','var')
    sclfactor=1;
end

[nparc,numClusters] = size(centroids);

if nparc > 400
    load('yeo7netlabelsLaus250.mat'); network7labels = network7labels(1:nparc);
else
    load('yeo7netlabelsLaus125.mat'); network7labels = network7labels(1:nparc);
end

binaryNetVectors = zeros(nparc,14);

for I = 1:7
    binaryNetVectors(:,I) = sclfactor*double(network7labels == I);
    binaryNetVectors(:,I+7) = -sclfactor*double(network7labels == I);
end

YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-'};

% calculate correlation between each binary state vector and each centroid

net7ED = zeros(numClusters,14);

for K = 1:numClusters
    for B = 1:14
        net7ED(K,B) = corr(centroids(:,K),binaryNetVectors(:,B));
    end
end

% get index of minimum and assign names

clusterNamesInit = cell(numClusters,1);
plusminus = true(numClusters,1);
for K = 1:numClusters
    ind = find(net7ED(K,:) == max(net7ED(K,:)));    % find maximum correlation
    clusterNamesInit{K} = YeoNetNames{ind};
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
