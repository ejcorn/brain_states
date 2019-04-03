function [G,Gc] = fcn_group_bins(A,dist,hemiid,nbins)
% fcn_distance_dependent_threshold      generate group matrix
%
%   G = fcn_distance_dependent_threshold(A,dist,hemiid,frac) generates a
%       group-representative structural connectivity matrix by preserving
%       within-/between-hemisphere connection length distributions.
%
%   Inputs:
%               A,      [node x node x subject] structural connectivity
%                       matrices.
%               dist,   [node x node] distance matrix
%               hemiid, indicator matrix for left (1) and right (2)
%                       hemispheres
%               nbins,  number of distance bins
%
%   Outputs:
%               G,      group matrix (binary) with distance-based consensus
%               Gc,     group matrix (binary) with traditional
%                       consistency-based thresholding.
%
%   Richard Betzel, Indiana University, 2018
%
% 

distbins = linspace(min(nonzeros(dist)),max(nonzeros(dist)),nbins + 1);
distbins(end) = distbins(end) + 1;

[n,~,nsub] = size(A); % number nodes (n) and subjects (nsub)
C = sum(A > 0,3);     % consistency
W = sum(A,3)./C;      % average weight
W(isnan(W)) = 0;      % remove nans
Grp = zeros(n,n,2);   % for storing inter/intra hemispheric connections (we do these separately)
Gc = Grp;
for j = 1:2
    if j == 1         % make inter- or intra-hemispheric edge mask
        d = +(hemiid == 1)*(hemiid' == 2);
        d = d | d';
    else
        d = +(hemiid == 1)*(hemiid' == 1) | +(hemiid == 2)*(hemiid' == 2);
        d = d | d';
    end
    m = dist.*d;
    D = nonzeros(bsxfun(@times,(A > 0),dist.*triu(d)));
    tgt = length(D)/nsub;
    G = zeros(n);
    for ibin = 1:nbins
        disp(ibin)  %EJC added
        mask = find(triu(m >= distbins(ibin) & m < distbins(ibin + 1),1));
        frac = round(tgt*sum(D >= distbins(ibin) & D < distbins(ibin + 1))/length(D));
        c = C(mask);
        [~,idx] = sort(c,'descend');
        G(mask(idx(1:frac))) = 1;
    end
    Grp(:,:,j) = G;
    I = find(triu(d,1));
    w = W(I);
    [~,idx] = sort(w,'descend');
    w = zeros(n);
    w(I(idx(1:nnz(G)))) = 1;
    Gc(:,:,j) = w;
end
G = sum(Grp,3); G = G + G';
Gc = sum(Gc,3); Gc = Gc + Gc';