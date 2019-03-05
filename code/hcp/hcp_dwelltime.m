addpaths;
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','hcpLR');
nparc = 462;
rTR = round(405); nTR = 405;

load(fullfile(savedir,['HCP_XHcentroids_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']),...
    'HCPcentroidsPNCorder','clusterNames','clusterNamesUp','clusterNamesDown','distanceMethod');
load(fullfile(savedir,['HCP_XHpartition_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']),...
    'HCPpartitionPNCorder','nsubjs','HCPsubjInd','HCPscanInd');

load(fullfile(savedir,['HCP_XHcentroids_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']));
load(fullfile(savedir,['HCP_XHpartition_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']));

scanlab = {'RestComb','nBackComb'};

for i = 1:numel(scanlab)
	TR = 0.72;     % HCP TR length
    HCPDwellTimeMean = zeros(nobs,numClusters);
    HCPDwellTimeMedian = zeros(nobs,numClusters);
    for N = 1:nobs
        [dt_mean,dt_median] = CALC_DWELL_TIME(HCPpartitionPNCorder(HCPsubjInd == N & HCPscanInd == (i-1)),numClusters);
        HCPDwellTimeMean(N,:) = dt_mean*TR;        % store dwell time in seconds
        HCPDwellTimeMedian(N,:) = dt_median*TR;        % store dwell time in seconds
    end	
	save(fullfile(savedir,[scanlab{i},'HCPDwellTime_k',num2str(numClusters),name_root,'.mat']),'HCPDwellTimeMean','HCPDwellTimeMedian');
end
