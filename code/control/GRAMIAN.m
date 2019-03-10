function [WcI] = GRAMIAN(A, T, normalize)
% Computes minimum control energy for state transition.
% A: System adjacency matrix: n x n
% x0: Initial state
% xf: Final state
% T: Control horizon
% 
% Outputs
% E: Minimum control energy
if ~exist('normalize','var')
	normalize = true;
end

% Normalize
if normalize
	A = (A / (max(eig(A)))) - eye(length(A));
	disp(['After normalization, max eigenvalue of A is ',num2str(max(eig(A)))])
end
% Gramian
Wc = integral(@(t)expm((A+A')*t), 0, T, 'ArrayValued', 1);
% Inverse
WcI = Wc^-1;