function Xo = L2MAGNITUDENORM(Xo)

% normalize Xo (NxM matrix that is a vector of M Nx1 states)
% by each states L2 magnitude

Xo_mag = sqrt(sum(Xo.^2,1));
Xo = Xo ./ Xo_mag;		% normalize all states by magnitude