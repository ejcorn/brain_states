addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['VolNormSC',num2str(lausanneScaleBOLD),'.mat']));

fname = ['/data/tesla-data/ecornblath/matlab/brainmapping2/LausanneNifti/ROIv_scale',num2str(lausanneScaleBOLD),'_dilated.nii.gz'];
lausnifti = load_nii(fname);
lauscoor = zeros(nparc,3);

for i = 1:nparc
    [xind,yind,zind] = ind2sub(size(lausnifti.img),find(ismember(lausnifti.img,i)));
    lauscoor(i,:) = mean([xind,yind,zind],1);
end

D = squareform(pdist(lauscoor,'Euclidean'));

SCvolnormNULL = cell(nobs,1);
for N = 1:nobs
    tic
    %A = randmio_und(SCvolnorm{N},1);
    nbins = 11; nrewire = 1e5; [~,A] = fcn_preserve_degseq_lengthdist(SCvolnorm{N},D,nbins,nrewire); 
    SCvolnormNULL{N} = A + A';    
    toc
    disp(['Subject = ',num2str(N)]);
end

save(['/data/tesla-data/ecornblath/matlab/control_fc/pipeline/data/VolNormSCRickNull',num2str(lausanneScaleBOLD),'.mat'],'SCvolnormNULL');