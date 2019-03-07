%zdim = 0; scan = 'C'; lausanneScaleBOLD = 250; extralabel = 'corrfinal';

name_root = ['Scan',scan,'Laus',num2str(lausanneScaleBOLD),'Z',num2str(zdim),extralabel];
datadir = [basedir,'data'];
masterdir = fullfile(basedir,'results',name_root);
demoLTN = readtable('/data/tesla-data/ecornblath/PNCLatestFreeze/cornblathPncRestNback/subject_demographics/n884_pnc_dtiNbackRest_subject_demographics.csv');
demoLTN = sortrows(demoLTN,'scanid','ascend');
nobs = length(demoLTN.scanid); demoLTN.Sex = demoLTN.sex - 1;
demoLTN.age_in_yrs = demoLTN.ageAtScan1 / 12;
rest_hm = demoLTN.restRelMeanRMSMotion;
nback_hm = demoLTN.nbackRelMeanRMSMotion;

brainvol = readtable([datadir,'/n1601_ctVol20170412.csv']);
brainvol = sortrows(brainvol,'scanid','ascend');
brainvol = brainvol(ismember(brainvol.scanid,demoLTN.scanid),:);
demoLTN.BrainSegVol = brainvol.mprage_antsCT_vol_TBV;

demo2 = readtable([datadir,'/n1601_dataRelease_MASTER_demogs.csv']);
demo2 = sortrows(demo2,'scanid','ascend');
demo2 = demo2(ismember(demo2.scanid,demoLTN.scanid),:);
demoLTN.handedness = demo2.handednessv2 - 1; % make handedness 0-1 binary

%% Load structure

if lausanneScaleBOLD == 250
    nparc = 462;
elseif lausanneScaleBOLD == 125
    nparc = 233;
end

dtifnames = cellstr([repmat(['/data/tesla-data/ecornblath/PNCLatestFreeze/cornblathPncRestNback/diffusion_data/volNormSC/Lausanne',num2str(lausanneScaleBOLD),'/'],[nobs 1]),num2str(demoLTN.scanid),...
    repmat(['_volNormStreamline_LausanneScale',num2str(lausanneScaleBOLD),'.mat'],[nobs 1])]);
    
% check if structure exists for every subject -- this removes 4 subjects in 462 node scale
StructureExist = logical(cellfun(@(x) exist(x,'file'),dtifnames));
demoLTN = demoLTN(StructureExist,:);
if sum(StructureExist) < nobs
    disp(['removing ',num2str(nobs-sum(StructureExist)),' subject(s) with missing structure'])
end
nobs = sum(StructureExist);
SCvolnorm = cell(nobs,1);
dtifnames = cellstr([repmat(['/data/tesla-data/ecornblath/PNCLatestFreeze/cornblathPncRestNback/diffusion_data/volNormSC/Lausanne',num2str(lausanneScaleBOLD),'/'],[nobs 1]),...
    num2str(demoLTN.scanid),...
    repmat(['_volNormStreamline_LausanneScale',num2str(lausanneScaleBOLD),'.mat'],[nobs 1])]);

% iterate through subjects, store structural connectivity matrices in a cell
for N = 1:nobs
    load(dtifnames{N});
    SCvolnorm{N} = volNorm_connectivity(1:nparc,1:nparc);  % remove brainstem  
end

disp('SC loaded');

%% Load BOLD data

restNumTRs = 120; nbackNumTRs = 225;

restTS = nan(restNumTRs,nparc,nobs);
nbackTS = nan(nbackNumTRs,nparc,nobs);

restfnames = [repmat(['data/tesla-data/ecornblath/PNCLatestFreeze/cornblathPncRestNback/rest_data/Lausanne',num2str(lausanneScaleBOLD),'/'],[nobs 1]),...
    num2str(demoLTN.scanid),...
    repmat(['_Lausanne',num2str(lausanneScaleBOLD),'_ts.1D'],[nobs 1])];
    restfnames = [repmat(['data/tesla-data/ecornblath/PNCLatestFreeze/cornblathPncRestNback/nback_data/Lausanne',num2str(lausanneScaleBOLD),'/'],[nobs 1]),...
    num2str(demoLTN.scanid),...
    repmat(['_Lausanne',num2str(lausanneScaleBOLD),'_ts.1D'],[nobs 1])];
nbackfnames = [repmat(['/data/jag/bassett-lab/Lausanne1601/nback/nbackNetwork_Lausanne',num2str(lausanneScaleBOLD),'/','Lausanne',num2str(lausanneScaleBOLD),'Timeseries/'],[nobs 1]),num2str(demoLTN.scanid),repmat(['_Lausanne',num2str(lausanneScaleBOLD),'_ts.1D'],[nobs 1])];

for N = 1:nobs
    disp(['Subject ',num2str(N)])
    % default: time series is NaNs
    tmpr = nan(restNumTRs,nparc);
    tmpn = nan(nbackNumTRs,nparc);
    
    % load time series and replace NaNs IF data exists
    if exist(restfnames{N},'file') && fseek(fopen(restfnames{N}), 1, 'bof') ~= -1
        tmpr = dlmread(restfnames{N});
    end
    if exist(nbackfnames{N},'file') && fseek(fopen(nbackfnames{N}), 1, 'bof') ~= -1
        tmpn = dlmread(nbackfnames{N});
    end
    
    % z-score according to specified input, if desired
    if zdim == 0
        restTS(:,:,N) = tmpr(:,1:nparc);
        nbackTS(:,:,N) = tmpn(:,1:nparc);
    elseif zdim == 1
        restTS(:,:,N) = zscore(tmpr(:,1:nparc),[],1);
        nbackTS(:,:,N) = zscore(tmpn(:,1:nparc),[],1);
    elseif zdim == 2
        restTS(:,:,N) = zscore(tmpr(:,1:nparc),[],2);
        nbackTS(:,:,N) = zscore(tmpn(:,1:nparc),[],2);
    end
end


disp('BOLD Data loaded');

% confirm there isn't any subject without missing data for either rest or n-back
rest_missingdata = ~logical(squeeze(sum(sum(isnan(restTS),1),2)));
nback_missingdata = ~logical(squeeze(sum(sum(isnan(nbackTS),1),2)));

% will be true only if subject has data for both rest AND n-back
BOLDMissingDataExclusionMask = and(rest_missingdata,nback_missingdata);

restTS = restTS(:,:,BOLDMissingDataExclusionMask);
nbackTS = nbackTS(:,:,BOLDMissingDataExclusionMask);
demoLTN = demoLTN(BOLDMissingDataExclusionMask,:);
rest_hm = rest_hm(BOLDMissingDataExclusionMask);
nback_hm = nback_hm(BOLDMissingDataExclusionMask);
nobs = sum(BOLDMissingDataExclusionMask);
SCvolnorm(find(~BOLDMissingDataExclusionMask)) = [];

%% concatenate

if strcmp(scan,'R')
    concTS = reshape(permute(restTS,[1 3 2]),nobs*restNumTRs,nparc);
    scanInd = zeros(nobs*restNumTRs,1);
    subjInd = cell2mat(arrayfun(@(a,r)repmat(a,1,r),1:nobs,repmat(restNumTRs,[1 nobs]),'uni',0));
elseif strcmp(scan,'N')
    concTS = reshape(permute(nbackTS,[1 3 2]),nobs*nbackNumTRs,nparc);
    scanInd = ones(nobs*nbackNumTRs,1);
    subjInd = cell2mat(arrayfun(@(a,r)repmat(a,1,r),1:nobs,repmat(nbackNumTRs,[1 nobs]),'uni',0));
elseif strcmp(scan,'C')
    concTS = [reshape(permute(restTS,[1 3 2]),nobs*restNumTRs,nparc);reshape(permute(nbackTS,[1 3 2]),nobs*nbackNumTRs,nparc)];
    scanInd = [zeros(nobs*restNumTRs,1);ones(nobs*nbackNumTRs,1)];
    subjInd = [cell2mat(arrayfun(@(a,r)repmat(a,1,r),1:nobs,repmat(restNumTRs,[1 nobs]),'uni',0)),cell2mat(arrayfun(@(a,r)repmat(a,1,r),1:nobs,repmat(nbackNumTRs,[1 nobs]),'uni',0))];
end

disp('Concatenation done');

%% Save data

cd(datadir);
save(['TimeSeriesIndicators',name_root,'.mat'],'scanInd','subjInd');
save(['VolNormSC',num2str(lausanneScaleBOLD),'.mat'],'SCvolnorm');
save(['Demographics',name_root,'.mat'],'demoLTN','nparc','nobs','lausanneScaleBOLD','scan','zdim','rest_hm','nback_hm','basedir','datadir','masterdir');
disp('.mat data saved');
csvwrite(['ConcTSCSV_',name_root,'.csv'],concTS);
writetable(demoLTN,['Demographics',name_root,'.csv']);
disp('.csv data saved');
disp('data saved');