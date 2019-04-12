function [MI,HX,HXY] = MUTUAL_INFO(X,Y)

% This script calculates the mutual information between integer sequences X and Y, i.e. I(X,Y) 
% 1. I(X,Y) = H(X) - H(X|Y)
% 2. H(X) = -sum_i[P(X = X_i)*log( P(X = X_i) )]	% base of log determines the units. base e = 'nats'
% 3. H(X|Y) = sum_i,j[P(X = X_i, Y = Y_j)*log( P(Y = Y_j) / P(X = X_i, Y = Y_j) )]
% 4. P(X,Y) = P(X=i & Y = j) for all pairs of i,j. 2D matrix summing to 1

X = double(X); Y = double(Y);
elem_X = unique(X);
elem_Y = unique(Y);
n_X = length(elem_X);
n_Y = length(elem_Y);
HX = zeros(n_X,1);
HXY = zeros(n_X,n_Y);
PX = zeros(n_X,1);
PY = zeros(n_Y,1);
PXY = zeros(n_X,n_Y);
for i = 1:n_X
	PX(i) = sum(X == elem_X(i)) / length(X);
	HX(i) = PX(i) * log(PX(i));
	for j = 1:n_Y
        PY(j) = sum(Y == elem_Y(j)) / length(Y);
		PXY(i,j) = sum(X == elem_X(i) & Y == elem_Y(j)) / length(X);
		HXY(i,j) = PXY(i,j) * log(PY(j)/PXY(i,j));
	end
end

MI = -sum(HX,'omitnan') - sum(sum(HXY,'omitnan'),'omitnan');
MI = MI / -sum(HX,'omitnan'); %Normalize to entropy of X