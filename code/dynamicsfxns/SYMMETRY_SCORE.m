function score = SYMMMETRY_SCORE(A)
	% sum(abs( off-diagonal of 0.5*(A-A')) / sum(abs( off-diagonal of A))
	% symm = 0, asymmetric = 1
	A = squeeze(A);
	n = size(A,1);
	score = sum(sum(abs(~eye(n) .* 0.5 * (A-A')))) / sum(sum(abs(~eye(n)*A)));
