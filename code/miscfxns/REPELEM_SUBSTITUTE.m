function Xr = REPELEM_SUBSTITUTE(X,N)

% repeat elements of row vector X n times

Xr = repmat(X',1,N)';
Xr=Xr(:)';