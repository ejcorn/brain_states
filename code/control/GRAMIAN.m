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
	A = (A / max(eigs(A,1)+1)) - eye(length(A));
end
% Gramian
Wc = integral(@(t)expm((A+A')*t), 0, T, 'ArrayValued', 1);
% Inverse
WcI = Wc^-1;