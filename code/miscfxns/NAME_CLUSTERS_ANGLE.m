function [clusterNames,reorderClusters,clusterNamesSort] = NAME_CLUSTERS_ANGLE(centroids,sclfactor)

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

% calculate E.D. from binary state vector to centroids

net7ED = zeros(numClusters,14);

for K = 1:numClusters
    for B = 1:14
        net7ED(K,B) = dot(centroids(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

% get index of minimum and assign names

clusterNamesInit = cell(numClusters,1);
plusminus = true(numClusters,1);
for K = 1:numClusters
    ind = find(net7ED(K,:) == max(net7ED(K,:)));    %cos(0) = 1 so need max not min (like for E.D.)
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
