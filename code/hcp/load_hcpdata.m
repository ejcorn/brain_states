% use this script on local machine to concatenate resting state and n-back HCP data into a data matrix concTS

basedir = '~/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states';
%% concatenate time series
rfilenames = dir('~/Dropbox/HCP_XH/rest/*LR.txt');
nfilenames = dir('~/Dropbox/HCP_XH/wm/*LR.txt');

nsubjs = 50; nparc = 462;
rTR = round(405); nTR = 405;
nparc = 462; nsubjs = length(rfilenames);
%
concTS = zeros((rTR+nTR)*nsubjs,nparc);

for subj = 1:nsubjs
    rTS = dlmread(fullfile(rfilenames(subj).folder,rfilenames(subj).name));
    st = 1+(subj-1)*rTR;
    concTS(st:(subj*rTR),:) = rTS(1:rTR,1:nparc);
    nTS = dlmread(fullfile(nfilenames(subj).folder,nfilenames(subj).name));
    st = (nsubjs*rTR)+ 1 + (subj-1)*nTR;
    concTS(st:(subj*nTR)+nsubjs*rTR,:) = nTS(1:nTR,1:nparc); 
end

save(fullfile(basedir,'data',['HCP_XH_concTS_R',num2str(rTR),'N',num2str(nTR),'.mat']),'concTS');