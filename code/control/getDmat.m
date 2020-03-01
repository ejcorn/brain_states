addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy');
mkdir(savedir);

%% get coordinates of spatial embedding

[nifti,sc_indices] = RETURN_NII_SUBCORT(name_root,lausanneScaleBOLD);

lauscoor = zeros(nparc,3);

for i = 1:nparc
    [xind,yind,zind] = ind2sub(size(nifti),find(ismember(nifti,i)));
    lauscoor(i,:) = mean([xind,yind,zind],1);
end

D = squareform(pdist(lauscoor,'Euclidean'));

save(fullfile(datadir,['Lausanne',num2str(lausanneScaleBOLD),'DistanceMatrix.mat']),'D');