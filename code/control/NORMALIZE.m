function Anorm = NORMALIZE(A)

% scale adjacency matrix in various ways

% normalize A to max eig +1 --> new max eig [0,1] - identity --> go to 0 over time
%Anorm = (A / (1 + max(eig(A))) ) - eye(length(A));

% normalize A to max eig --> new max eig = 1 - identity --> go to single dominant eigenvector over time
Anorm = (A / (max(eig(A))) ) - eye(length(A));
