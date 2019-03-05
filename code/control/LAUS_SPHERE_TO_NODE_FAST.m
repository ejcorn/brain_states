function [data] = LAUS_SPHERE_TO_NODE_FAST(Llabels,Rlabels,nparc,annot)

% this function takes sphere coordinates and converts it back to node-level data
Lv = annot.Lv; LL = annot.LL; Lct = annot.Lct; Rv = annot.Rv; RL = annot.RL; Rct = annot.Rct;

if ~exist('nparc','var')
	nparc = 462;	% lausanne 250
end

if nparc == 233
    lausanneScale = 125; idx = 3; mid = 115;
elseif nparc == 462
    lausanneScale = 250; idx = 4; mid = 230;
end

lausnames = annot.roinames{idx};

%%

Ls = size(Lct.struct_names,1); Rs = size(Rct.struct_names,1);

for i = 1:Ls
    Lct.struct_names{i} = ['lh_',Lct.struct_names{i}];
end

for i = 1:Rs
    Rct.struct_names{i} = ['rh_',Rct.struct_names{i}];
end

[~,lind] = ismember(Lct.struct_names,lausnames);
[~,rind] = ismember(Rct.struct_names,lausnames);

% lind and rind will contain the ROI indices in lausnames, ordered as in
% ?ct.struct_names. i.e. the ROI indices you work with outside of
% annotation files, in the order they are found in the color table of the
% annotation file
% this is the key step in converting data ordered the usual lausanne way
% to surface coordinates
% this is not an eloquent way of explaining that fact
% Below I will call the usual lausanne ROI indices 'MatrixIndices'
% because that's how they're ordered in structural adjmats
% I will call annotation index 'AnnotationIndices'

%% now assign all the values in labels to nodes

data = zeros(nparc,size(Rlabels,2));	% nparc by permno in SpinPermuFS.m essentially. if Rlabels is vector then it's just (nparc,1)

L_labelKey = Lct.table(:,5);
for i = 1:length(lind)
	MatrixIndex = lind(i);
	if MatrixIndex ~= 0		% don't use AnnotationIndices that refer to regions that don't exist, like medial wall
		data(MatrixIndex,:) = nanmean(Llabels(LL == L_labelKey(i),:),1);
	end
end

R_labelKey = Rct.table(:,5);
for i = 1:length(rind)
	MatrixIndex = rind(i);
	if MatrixIndex ~= 0		% don't use AnnotationIndices that refer to regions that don't exist, like medial wall
		data(MatrixIndex,:) = nanmean(Rlabels(RL == R_labelKey(i),:),1);
	end
end

