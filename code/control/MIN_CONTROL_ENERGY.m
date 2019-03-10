function [E] = MIN_CONTROL_ENERGY(A, WcI, x0, xf, T,normalize)
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
% State transition to achieve
Phi = expm(A*T)*x0 - xf;
% Energy
E = sum((WcI*Phi).*Phi);
end