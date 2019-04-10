function Anorm = NORMALIZE(A,c)

% scale adjacency matrix in various ways

% normalize A to max eig + c where c >0 --> new max eig [0,1] - identity --> go to 0 over time
% normalize A to max eig --> new max eig = 1 - identity --> go to single dominant eigenvector over time

if ~exist('c','var')
	disp('automatically setting c=0 so that A is normalized to a max eigenval of 0.')
	c=0;
end
Anorm = (A / (c+ max(eig(A))) ) - eye(length(A));
