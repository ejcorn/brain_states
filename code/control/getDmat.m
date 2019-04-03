addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','control_energy');
mkdir(savedir);

%% get coordinates of spatial embedding

fname = ['data/ROIv_scale',num2str(lausanneScaleBOLD),'_dilated.nii.gz'];
lausnifti = load_nii(fname);
lauscoor = zeros(nparc,3);

for i = 1:nparc
    [xind,yind,zind] = ind2sub(size(lausnifti.img),find(ismember(lausnifti.img,i)));
    lauscoor(i,:) = mean([xind,yind,zind],1);
end

D = squareform(pdist(lauscoor,'Euclidean'));

save(fullfile(datadir,['Lausanne',num2str(lausanneScaleBOLD),'DistanceMatrix.mat']),'D');