addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load([masterdir,'/clusterAssignments/repkmeansPartitions',distanceMethod,name_root,'.mat']);
savedir = fullfile(masterdir,'analyses','choosing_k');
mkdir(savedir);

minNumClusters = 2; maxNumClusters = 11; clusterRange = minNumClusters:maxNumClusters; numK = 1 + maxNumClusters-minNumClusters;
k_offset = minNumClusters - 1;

for numClusters = minNumClusters:maxNumClusters
	f=figure;
	histogram(combSumD{numClusters-k_offset});
	xlabel('Cluster Solution Error');
	ylabel('# of Solutions');
	title(['Cluster Error Distribution, k=',num2str(numClusters)]);
	set(gca,'FontSize',8);
	f.PaperUnits = 'inches';
	f.PaperPosition = [0 0 3 3];
	f.PaperSize = [3 3];
	xlim(get(gca,'XLim'));		% don't change x limits after resizing figure, keep automatic ones
	ylim([0 length(combSumD{numClusters-k_offset})]);	% 0 to number of initializations
	saveas(f,fullfile(savedir,['ClusterErrorDistribution_k',num2str(numClusters),'.pdf']));
end