function [WcI] = GRAMIAN_FAST(A, T)
% Computes inverse controllability Gramian for continuous systems
% Input
% A: System adjacency matrix, MUST BE SYMMETRIC: n x n
% T: Time horizon: Should be < 1e6
% 
% Output:
% WcI: Inverse of Controllability Gramian

% Gramian
[V, U] = eig(A);
u = diag(U);

% Find vanishingly small eigenvalues
sInd = find(abs(u) < 1e-10);        % Small Magnitude Eigenvalue Indices
uInd = setdiff(1:size(A,1), sInd);  % Larger MagnitudeEigenvalue Indices

% Invert Small Values
WcI = V(:,uInd)*diag((2*u(uInd))./(exp(2*u(uInd)*T)-1))*V(:,uInd)';

% Add Large Values
WcI = WcI + V(:,sInd)*V(:,sInd)'/T;