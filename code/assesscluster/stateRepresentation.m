load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));

scanlab = {'RestComb','nBackComb'};
minNumClusters = 2; maxNumClusters = 11; 
clusterRange = minNumClusters:maxNumClusters; 
numK = 1 + maxNumClusters-minNumClusters;

restAbsentStateCount = zeros(1,numK);
nbackAbsentStateCount = zeros(1,numK);
dwelltime = cell(numK,2);
for C = 1:numK
	numClusters = clusterRange(C);
	load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
	kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;
	
	for S = 1:numel(scanlab)
		dwelltime{C,S} = zeros(nobs,numClusters);
		tmpAssignments = kClusterAssignments(scanInd == (S-1));
		tmpSubjInd = subjInd(scanInd == (S-1));
		for N = 1:nobs
			subjStates = tmpAssignments(tmpSubjInd == N);
			for K = 1:numClusters
				dwelltime{C,S}(N,K) = sum(subjStates == K);
			end
		end
	end
	restAbsentStateCount(C) = 100* sum(sum(dwelltime{C,1} == 0,2) > 0) / nobs;	% percent of subjects with >0 missing states
	nbackAbsentStateCount(C) = 100* sum(sum(dwelltime{C,2} == 0,2) > 0) / nobs;
end

RNcolors = {'#005C9F','#FF8400'};
RNcolors = hex2rgb(RNcolors);
f = figure;
bar([restAbsentStateCount;nbackAbsentStateCount]','FaceAlpha',.5);
xticks(1:numK);xticklabels(minNumClusters:maxNumClusters);
ylabel('% Subjects'); xlabel('k','FontWeight','bold'); title('Absent States');
legend({'Rest','n-back'},'Location','northwest')
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

savedir = fullfile(masterdir,'analyses','choosing_k');
mkdir(savedir);

saveas(f,fullfile(savedir,['StateRepresentation',name_root,'.pdf']));
CombinedStateRepresentation = [restAbsentStateCount;nbackAbsentStateCount]';
save(fullfile(savedir,['FigS1b__StateRepresentation',name_root,'.mat']),'CombinedStateRepresentation');
