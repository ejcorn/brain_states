function A_DLWnull = DLW_NULL(A,D)

% wrapper function to apply null model from Betzel et al. 2018 Specificity and robustness of long-distance connections
% A: connectivity matrix
% D: euclidean distance between regions
nbins = 11; nrewire = 1e5; 	% parameters for spatial null model
[~,A_DLWnull] = fcn_preserve_degseq_lengthdist(A,D,nbins,nrewire); 
A_DLWnull = A_DLWnull + A_DLWnull';