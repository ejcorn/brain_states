function [VarianceExplained,WithinClusterVariance,BetweenClusterVariance] = VAREXPLAINED(X,partition,Centroids,numClusters,distanceMethod)

% X: N-by-P data matrix
% Centroids: P-by-numClusters centroids
% partition: vector of cluster indices obtained by clustering X
% numClusters: number of clusters
% distanceMethod: function for computing distance in k-means

if ~exist('distanceMethod','var')
    distanceMethod = 'correlation'; % default distance function is the one we used in the paper (correlation)
end

if distanceMethod == 'correlation'
    dfun = @(x,y) (1-corr(x,y)).^2;
elseif distanceMethod == 'sqeuclidean'
    dfun = @(x,y) sum((x-y).^2,1);
end

N = length(partition);
overallDataCentroid = zeros(size(Centroids)); % overall data centroid is average of cluster centroids weighted by number of observations in each cluster
for k = 1:numClusters
    overallDataCentroid(:,k) = sum(partition == k)*Centroids(:,k); % weighted sum of each cluster centroid
end
overallDataCentroid = (1/N)*sum(overallDataCentroid,2); % divide by number of observations to get weighted average

WithinClusterVariance = zeros(numClusters,1);
BetweenClusterVariance = zeros(numClusters,1); % average distance of centroids to overall center of data, as proxy for average distance of all points to center of space
for k = 1:numClusters
    WithinClusterVariance(k) = sum(dfun(Centroids(:,k),X(partition == k,:)')); % distance between each point and its cluster center
    BetweenClusterVariance(k) = sum(partition == k)*sum(dfun(Centroids(:,k),overallDataCentroid)); % sq. distance between each cluster centroid and overall data centroid, weighted by number of points in each centroid
end
WithinClusterVariance = sum(WithinClusterVariance) / N; % compute average distance of any sample to its cluster centroid
BetweenClusterVariance = sum(BetweenClusterVariance) / N; % compute weighted average of distance between each cluster centroid and overall data centroid
VarianceExplained = BetweenClusterVariance / (WithinClusterVariance + BetweenClusterVariance);