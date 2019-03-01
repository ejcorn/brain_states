addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['VolNormSC',num2str(lausanneScaleBOLD),'.mat']));

SCvolnormNULL = cell(nobs,1);
for N = 1:nobs
    tic
    SCvolnormNULL{N} = randmio_und(SCvolnorm{N},1);;    
    toc
    disp(['Subject = ',num2str(N)]);
end

save(fullfile(datadir,['VolNormSCBCTNull',num2str(lausanneScaleBOLD),'.mat']),'SCvolnormNULL');