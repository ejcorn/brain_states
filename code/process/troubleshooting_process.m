% this script was written to compare the data set I used for the submitted manuscript with the data set I used for the revised manuscript

addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));

restfnames = cellstr([repmat(['/data/tesla-data/ecornblath/PNCLatestFreeze/cornblathPncRestNback/rest_data/Lausanne',num2str(lausanneScaleBOLD),'/'],[nobs 1]),...
    num2str(demoLTN.scanid),... % only load files in demoLTN
    repmat(['_Lausanne',num2str(lausanneScaleBOLD),'_ts.1D'],[nobs 1])]);

nbackfnames = cellstr([repmat(['/data/tesla-data/ecornblath/PNCLatestFreeze/cornblathPncRestNback/nback_data/Lausanne',num2str(lausanneScaleBOLD),'/'],[nobs 1]),...
    num2str(demoLTN.scanid),... % only load files in demoLTN
    repmat(['_Lausanne',num2str(lausanneScaleBOLD),'_ts.1D'],[nobs 1])]);

old_restfnames = [repmat(['/data/jag/bassett-lab/Lausanne1601/rest/restNetwork_Lausanne',num2str(lausanneScaleBOLD),'/','Lausanne',num2str(lausanneScaleBOLD),'Timeseries/'],[nobs 1]),num2str(demoLTN.scanid),repmat(['_Lausanne',num2str(lausanneScaleBOLD),'_ts.1D'],[nobs 1])];
old_nbackfnames = [repmat(['/data/jag/bassett-lab/Lausanne1601/nback/nbackNetwork_Lausanne',num2str(lausanneScaleBOLD),'/','Lausanne',num2str(lausanneScaleBOLD),'Timeseries/'],[nobs 1]),num2str(demoLTN.scanid),repmat(['_Lausanne',num2str(lausanneScaleBOLD),'_ts.1D'],[nobs 1])];

num_subjs = nobs;
id_rest = zeros(num_subjs,1);
id_nback = zeros(num_subjs,1);
r_rest = zeros(num_subjs,1);
r_nback = zeros(num_subjs,1);
for N = 1:num_subjs
	disp(['Subject ',num2str(N)])
	new_ts = dlmread(restfnames{N});
	old_ts = dlmread(old_restfnames(N,:));
	id_rest(N) = isequal(new_ts,old_ts);	% are resting state time series identical
	r_rest(N) = corr(reshape(new_ts,[],1),reshape(old_ts,[],1));
	new_ts = dlmread(nbackfnames{N});
	old_ts = dlmread(old_nbackfnames(N,:));
	id_nback(N) = isequal(new_ts,old_ts); % are n-back time series identical
	r_nback(N) = corr(reshape(new_ts,[],1),reshape(old_ts,[],1));

end

alt_subjs = find(id_rest == 0);	% get subjects with differences in rest
num_subjs = length(alt_subjs);
r_rest = zeros(num_subjs,nparc);
for N = 1:num_subjs
	disp(['Subject ',num2str(N)])
	subj = alt_subjs(N);
	new_ts = dlmread(restfnames{subj});
	old_ts = dlmread(old_restfnames(subj,:));
	for R = 1:nparc		% look at time series by region
		r_rest(N,R) = corr(new_ts(:,R),old_ts(:,R));
	end
end