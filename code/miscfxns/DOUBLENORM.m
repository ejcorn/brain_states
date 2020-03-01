function [Xo,Xf] = DOUBLENORM(Xo,Xf)

% normalize Xo and Xf (NxM matrices that are vectors of M Nx1 states)
% by their magnitude, and then by the distance between pairs of columns in Xo and Xf

Xo_mag = sqrt(sum(Xo.^2,1));
Xo = Xo ./ Xo_mag;		% normalize all states by magnitude
Xf_mag = sqrt(sum(Xf.^2,1));
Xf = Xf ./ Xf_mag;		% normalize all states by magnitude

XoXf_mag = sqrt(sum((Xo-Xf).^2,1)); % magnitude of interstate distance
Xo = Xo ./ XoXf_mag; % normalize by inter-state distance
Xf = Xf ./ XoXf_mag; % normalize by inter-state distance
