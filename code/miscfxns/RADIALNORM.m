function [x, z] = RADIALNORM(x,y)

% for pairs of columns in x and y, normalize y so that it is unit
% distance from x

% this is accomplished by finding the point z, which is equal to x + the
% unit vector in the direction of x-y, i.e. :
% z = x + (x-y) / sqrt(sum((x-y).^2))

unit_v = (y-x) ./ sqrt(sum((x-y).^2));
z = x + unit_v;

% check:
% disp(sum((z-x).^2))
% should be all 1's