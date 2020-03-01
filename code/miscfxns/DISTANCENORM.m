function [Xo,Xf,XoXf_mag] = DISTANCENORM(Xo,Xf)

% normalize Xo and Xf (NxM matrices that are vectors of M Nx1 states)
% by the distance between pairs of columns in Xo and Xf

XoXf_mag = sqrt(sum((Xo-Xf).^2,1)); % magnitude of interstate distance
Xo = Xo ./ XoXf_mag; % normalize by inter-state distance
Xf = Xf ./ XoXf_mag; % normalize by inter-state distance
