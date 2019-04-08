addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));

for N = 1:nobs
	disp(['Subject ', num2str(N)]);
    % rest independent phase randomization
    concTS(subjInd' == N & scanInd == 0,:) = linsurr_ind(concTS(subjInd' == N & scanInd == 0,:));
    % nback independent phase randomization
    concTS(subjInd' == N & scanInd == 1,:) = linsurr_ind(concTS(subjInd' == N & scanInd == 1,:));    
end
iprTS = concTS;

save(fullfile(datadir,['iprtimeseries',name_root,'.mat']),'iprTS');
clear concTS

randTS = randn(size(iprTS));
save(fullfile(datadir,['randtimeseries',name_root,'.mat']),'randTS');
