function [runs,n_runs] = GET_BLOCK_STATE_TIME_SERIES(partition,indicator)

% partition: Vector of integers corresponding to state time series
% indicator: binary vector indicating which states belong to which blocks

% given a partition with possible entries 1:k
% return state time series for each block specified in indicator
indicator = int8(indicator);
indicator = reshape(indicator,1,[]);	% make partition row vector
run_start_stop = find([1,diff(indicator),1]);	% get indices of start and stops of each run
n_runs = length(run_start_stop)-1;	% total number of runs
runs = cell(1);
% iterate through runs, store each run for indicator == 1 in a cell array
for j = 1:n_runs	
	run_start = run_start_stop(j);
	run_stop = run_start_stop(j+1)-1;
	run_idx = sum(~cellfun(@isempty, runs)) + 1;
	if sum(indicator(run_start:run_stop)) ~= 0		% if the run corresponds to a block selected by the binary indicator
		runs{run_idx} = partition(run_start:run_stop);		% store that run
	end
end

n_runs = numel(runs);		% how many runs in the block