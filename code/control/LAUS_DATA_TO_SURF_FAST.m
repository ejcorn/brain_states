function [Rlabels,Llabels] = LAUS_DATA_TO_SURF_FAST(data,annot)

% this function takes lausanne node data and maps to fsaverage surface
Lv = annot.Lv; LL = annot.LL; Lct = annot.Lct; Rv = annot.Rv; RL = annot.RL; Rct = annot.Rct;
Lvertices = annot.Lvertices; Lfaces = annot.Lfaces; Rvertices = annot.Rvertices; Rfaces = annot.Rfaces;

nparc = size(data,1);
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

%convert names to numbers
%make true the indices that match your betas
%replace the trues with the beta values

%% 
Lids = false(1,Ls); Rids = false(1,Rs); indices = 1:nparc;
Lids(find(ismember(lind,indices))) = true; Rids(find(ismember(rind,indices))) = true; 
%Lids is a logical the size of Lct.table (contains node labels, reference
%for LL, which labels vertices). Lids tells you which nodes in indices are
%on the left brain.
k = 1; Linds = []; Lbetas = [];
for i = 1:nparc
    if ismember(i,lind)
        Linds(k,1) = Lct.table(find(ismember(lind,i)),5);
        Lbetas(k,1) = data(i);
        k = k + 1;
    end
end

k = 1; Rinds = []; Rbetas = [];
for i = 1:nparc
    if ismember(i,rind)
        Rinds(k,1) = Rct.table(find(ismember(rind,i)),5);
        Rbetas(k,1) = data(i);
        k = k + 1;
    end
end

%% get faces for all regions not contained in indices in grey

%Lvertices(ismember(LL,Linds)); Rvertices(ismember(RL,Rinds)); %all vertices corresponding to areas to remove
Lremvertind = find(ismember(LL,Linds)); Rremvertind = find(ismember(RL,Rinds)); %indices of all vertices to remove
for i = 1:3
    Lx(:,i) = ismember(Lfaces(:,i),Lremvertind);
    Rx(:,i) = ismember(Rfaces(:,i),Rremvertind);
end
% every entry of Lx/Rx that is a 1 denotes an index in Lfaces/Rfaces
% referencing a vertex for a region that needs to be deleted
 %rows of Lx/Rx with any 1s in them should be deleted because they belong to
 %remvertind (refer to vertices of ROIs in indices/betas)
 
Lfacesnull = Lfaces; Lfacesnull(any(Lx,2),:) = [];
Rfacesnull = Rfaces; Rfacesnull(any(Rx,2),:) = [];

%% get faces for all regions contained in indices on a color scale

Lremvertind = find(~ismember(LL,Linds)); Rremvertind = find(~ismember(RL,Rinds)); %indices of all vertices to remove
%any vertex not corresponding to a member of Linds must be removed if you
%only want to plot regions in indices
for i = 1:3
    Lx(:,i) = ismember(Lfaces(:,i),Lremvertind);
    Rx(:,i) = ismember(Rfaces(:,i),Rremvertind);
end
Lfaces(any(Lx,2),:) = [];
Rfaces(any(Rx,2),:) = [];

Llabels = zeros(size(LL)); Rlabels = zeros(size(RL));


for i=1:length(Linds)
    Llabels(LL==Linds(i)) = Lbetas(i);
end

for i=1:length(Rinds)
    Rlabels(RL==Rinds(i)) = Rbetas(i);
end
