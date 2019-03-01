a = clock;
rng(a(6));
addpaths;
minNumClusters = 2; maxNumClusters = 18; clusterRange = minNumClusters:maxNumClusters; numK = 1 + maxNumClusters-minNumClusters;
k_offset = minNumClusters - 1;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load([masterdir,'/clusterAssignments/repkmeansPartitions',distanceMethod,name_root,'.mat']);

savedir = [masterdir,'/clusterAssignments'];
%% calculate zrand
nCombos = nreps*(nreps-1)/2;
combos = nchoosek(1:nreps,2);
combos1 = combos(:,1); combos2 = combos(:,2);
%zRandPartitions = cell(numK,1);
zRandK = zeros(numK,2);
%parpool(maxNumCompThreads)
K = numClusters;
tmp = zeros(nCombos,1);
for i = 1:nCombos
    disp(i)	
    part1 = combos1(i); part2 = combos2(i);
    tmp(i) = zrand(double(combPartitions{K - k_offset}(:,part1)),double(combPartitions{K - k_offset}(:,part2)));
end
zRandPartitions = tmp;
zRandK(K,1) = mean(zRandPartitions); zRandK(K,2) = std(zRandPartitions);
disp(['zrand ',num2str(K)]);

%% plot
cd(savedir);
save(['zrand_k',num2str(numClusters),name_root,'.mat'],'zRandK');
