function [nifti,sc_indices,hemi] = RETURN_NII_SUBCORT(atlasName,atlasScale)

% INPUTS:
% atlasName: string containig atlas name
% atlasScale: number of parcels or scale identifier
%
% OUTPUTS:
% nifti: nifti file for atlas and scale
% sc_indices: indices of subcortical ROIs
% hemi: indices of hemisphere. 1 for left hemi, 2 for right hemi

if contains(atlasName,'Laus')
    nifti = load_nii(['data/nifti/ROIv_scale',num2str(atlasScale),'.nii.gz']);
    nifti = nifti.img;

    if atlasScale == 125
        sc_indices = [109:115,227:233];
        hemi = [2*ones(115,1);ones(118,1)];
    elseif atlasScale == 250
        sc_indices = [224:230,456:462];
        hemi = [2*ones(230,1);ones(232,1)];
    end
elseif contains(atlasName,'Brainnetome')
    nifti = load_nii(['data/nifti/BN_Atlas_',num2str(atlasScale),'_2mm.nii.gz']);
    nifti = nifti.img;

    if atlasScale == 246
        sc_indices = [211:246];
        hemi = repmat([1 2]',[atlasScale/2 1]);
    end
elseif contains(atlasName,'Schaefer')
    nifti = load_nii(['data/nifti/Schaefer2018_',num2str(atlasScale),'Parcels_7Networks_order_FSLMNI152_2mm.nii.gz']);
    nifti = nifti.img;
    sc_indices = []; % no subcortex
    hemi = [ones(atlasScale/2,1);2*ones(atlasScale/2,1)];
end
