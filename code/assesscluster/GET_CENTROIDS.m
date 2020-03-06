function kClusterCentroids = GET_CENTROIDS(X,partition,k)

% INPUTS:
% X: N observations by p features data matrix
% partition: length N integer partition vector 
% k: number of clusters in partition
%
% OUTPUTS:
% kClusterCentroids: pxk matrix whose columns contain the average value of points in each cluster
p = size(X,2);
kClusterCentroids = zeros(p,k);

for j = 1:k
    kClusterCentroids(:,j) = mean(X(partition == j,:),1)';
end