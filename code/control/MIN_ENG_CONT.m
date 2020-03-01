function [ x, u, n_err ] = MIN_ENG_CONT(A, T, B, x0, xf, nor)
% Computes minimum control energy for state transition.
% A: System adjacency matrix:         N x N
% B: Control input matrix:            N x k
% x0: Initial state:                  N x 1
% xf: Final state:                    N x 1
% T: Control horizon                  1 x 1
% nor: normalization                  boolean
% 
% Outputs
% x: State Trajectory
% u: Control Input

% Normalize
if nor == 1
    A = A/(eigs(A,1)+1) - eye(size(A));
end

% System Size
n = size(A,1);

% Compute Matrix Exponential
AT = [A              -.5*(B*B');...
      zeros(size(A)) -A'];
E = expm(AT*T);

% Compute Costate Initial Condition
E12 = E(1:n,[1:n]+n);
E11 = E(1:n,1:n);
p0 = pinv(E12)*(xf - E11*x0);

% Compute Costate Initial Condition Error Induced by Inverse
n_err = norm(E12*p0 - (xf - E11*x0));

% Prepare Simulation
nStep=1000;
t = linspace(0,T,nStep+1);

v0 = [x0; p0];                      % Initial Condition
v = zeros(2*n,length(t));           % Trajectory
Et = expm(AT*T/(length(t)-1));
v(:,1) = v0;

% Simulate State and Costate Trajectories
for i = 2:length(t)
    v(:,i) = Et*v(:,i-1);
end
x = v(1:n,:);
u = -0.5*B'*v([1:n]+n,:);

disp([n_err, norm(x(:,end)-xf)]);

% transpose to be similar to opt_eng_cont
u = u';
x = x';

end

