function [dt,run_lengths] = CALC_DWELL_TIME(partition,k)

% given a partition with possible entries 1:k
% find the length of runs of each cluster

partition = reshape(partition,1,[]);	% make partition row vector
run_start_stop = find([1,diff(partition),1]);	% get indices of start and stops of each run
n_runs = length(run_start_stop)-1;	% total number of runs
[temp_runs{1:k,1}] = deal(cell(1));
% iterate through runs, store in a cell array of cells for each cluster
% each cluster is a cell, then cells within those cells hold runs
for j = 1:n_runs
	run_start = run_start_stop(j);
	run_stop = run_start_stop(j+1)-1;
	run_j = partition(run_start:run_stop);
	% 
	run_k = unique(run_j);	% which cluster does the run belong to
    % find index to replace initial empty entry or add onto existing entries 
    run_k_idx = sum(~cellfun(@isempty, temp_runs{run_k})) + 1;
	temp_runs{run_k}{run_k_idx} = run_j;
end

dt = zeros(k,1);
run_lengths = cell(k,1);
for k_i = 1:k
    % store mean length of runs for each cluster
    run_lengths{k_i} = cellfun(@length,temp_runs{k_i});
    dt(k_i) = mean(run_lengths{k_i});
end
